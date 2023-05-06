%% Constellation macro-orbit determination 
% Date: 01/30/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of an Orbit object

classdef Orbit
    % Fundamental definition 
    properties 
        Label = "MyOrbit";
        mu                      % Central body gravitational parameter

        J2 = 0;                 % First zonal harmonic of the central body 
        Re = 0;                 % Reference radius of the central body for its spherical Legendre gravitational expansion
    
        ElementType             % Type of orbital elements in use
        ElementSet              % Parametrizing orbital element set
        Sigma                   % Covariance matrix of the estimated element set

        InitialEpoch            % Initial epoch 
        CurrentEpoch            % Current propagation epoch 
        PropagatedEpoch;        % End time propagation epoch
        FinalEpoch              % Final orbit epoch
        TimeStep                % Integration time step

        StateEvolution = [];    % Evolution of the orbital element set
        Dynamics;               % Dynamical model to be used

        % Gravitational parameters of the Sun and Moon
        mus = 1.32712440018e20;
        mul = 4.9048695e12;
    end

    % Private properties 
    properties (Hidden = true)
        m = 6;                  % State space dimension

        Tspan                   % Integration time span

        IntegrationOptions;     % Integration toleranc
        Trajectory;             % Cartesian trajectory

        Lc;                     % Characteristic length 
        Tc;                     % Characteristic time 
        Normalized = false;     % Flag to account for normalized coordinates
    end

    % Public methods 
    methods
        % General constructor
        function [obj] = Orbit(mu, myElementType, myElementSet, myInitialEpoch)
            % User specified values
            obj.ElementSet = myElementSet;
            obj.InitialEpoch = myInitialEpoch;
            obj.CurrentEpoch = obj.InitialEpoch;
            obj.PropagatedEpoch = obj.InitialEpoch;
            obj.FinalEpoch = obj.InitialEpoch;

            obj.mu = mu;

            % Dafault values 
            obj.IntegrationOptions = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      % Integration tolerances

            % Sanity checks 
            if (length(myElementSet) < obj.m)
                error('The state vector dimension is not valid.');
            end

            if (size(myElementSet,2) <= 1)
                myElementSet = myElementSet.';
            end

            % Assignment
            switch (myElementType)
                case 'ECI'
                    obj.ElementType = myElementType;
                case 'POL'
                    obj.ElementType = myElementType;
                case 'COE'
                    obj.ElementType = myElementType;
                case 'MOE'
                    obj.ElementType = myElementType;
                case 'KS'
                    obj.ElementType = myElementType;
                otherwise
                    error('No valid stated vector was introduced.')
            end

            try 
                switch (myElementType)
                    case 'ECI'
                        obj.ElementSet = myElementSet(1,1:obj.m);
                    case 'POL'
                        obj.ElementType = myElementType;
                    case 'COE'
                        obj.ElementSet = myElementSet(1,1:obj.m);
                        obj.ElementSet = [obj.ElementSet obj.ElementSet(1)*(1-obj.ElementSet(2)^2)];
                    case 'MOE'
                        obj.ElementSet = myElementSet(1,1:obj.m);
                    case 'KS'
                        obj.ElementSet = myElementSet(1,1:obj.m+2);
                    otherwise
                        error('No valid stated vector was introduced.')
                end
            catch ME
                rethrow(ME);
            end

            obj.StateEvolution(1,:) = [0, obj.ElementSet];
        end

        % Change the orbit label
        function [obj] = ChangeLabel(obj, myName)
            obj.Name = myName;
        end

        % Add final epoch 
        function [obj] = SetFinalEpoch(obj, myFinalEpoch)
            if (myFinalEpoch < obj.InitialEpoch)
                error('Specified final epoch is not greater than the initial epoch');
            else
                obj.FinalEpoch = myFinalEpoch;
            end
        end


        % Specify current epoch 
        function [obj] = SetCurrentEpoch(obj, myCurrentEpoch)
            if (myCurrentEpoch < obj.CurrentEpoch)
                error('New current epoch is not greater than the actual current epoch');
            elseif (myCurrentEpoch <= obj.FinalEpoch)
                obj.CurrentEpoch = myCurrentEpoch;
            elseif (myCurrentEpoch > obj.FinalEpoch)
                obj.FinalEpoch = myCurrentEpoch;
                obj.CurrentEpoch = myCurrentEpoch;
            end
        end

        % Add covariance 
        function [obj] = AddCovariance(mySigma)
            obj.Sigma = mySigma; 

            % Sanity check 
            try chol(mySigma)
            catch ME
                error('Covariance matrix is not positive definite.');
            end
        end

        % Propagation tolerances 
        function [obj] = AddPropagator(obj, myModel, myTimeStep)
            % Add the timestep
            if (myTimeStep < (obj.CurrentEpoch-obj.InitialEpoch)/2)
                obj.TimeStep = (obj.CurrentEpoch-obj.InitialEpoch)/2;
            else
                obj.TimeStep = myTimeStep;
            end

            % Select the model
            obj.Dynamics = obj.AddModel(myModel);
        end

        % Add integration tolerances
        function [obj] = AddIntegrationTolerances(myIntegrationOptions)
            obj.IntegrationOptions = myIntegrationOptions;
        end

        % Define J2 problem 
        function [obj] = DefineJ2Problem(obj, myJ2, myRe)
            obj.J2 = myJ2; 
            
            if (obj.Normalized)
                obj.Re = myRe/obj.Lc;
            else
                obj.Re = myRe;
            end
        end
        
        % Propagate to current epoch 
        function [obj] = Propagate(obj)
            % Check time span limits 
            if (obj.PropagatedEpoch < obj.CurrentEpoch)
                aux_tspan = 0:obj.TimeStep:(obj.CurrentEpoch-obj.PropagatedEpoch) * 24 * 3600;

                if (obj.Normalized)
                    aux_tspan = aux_tspan / obj.Tc;
                end
    
                % Perform the propagation 
                AuxEvolution = obj.Dynamics(aux_tspan); 

                % Update states
                obj.StateEvolution = [obj.StateEvolution; [obj.StateEvolution(end,1)+aux_tspan(2:end).' AuxEvolution(2:end,:)]];
                obj.ElementSet = AuxEvolution(end,:);

                obj.Tspan = [obj.Tspan obj.PropagatedEpoch:obj.TimeStep / (24 * 3600):obj.CurrentEpoch];
    
                % Update the last propagated epoch 
                obj.PropagatedEpoch = obj.CurrentEpoch;
            end
        end

        % Change the state evolution
        function [obj] = ChangeStateFormat(obj, myElementType)
            switch (myElementType)
                case obj.ElementType
                case 'ECI'
                    obj = State2Cartesian(obj);
                case 'POL'
                    obj = State2POL(obj);
                case 'COE'
                    obj = State2COE(obj);
                case 'MOE'
                    obj = State2MOE(obj);
                case 'KS'
                    obj = State2KS(obj);
                otherwise
                    error('No valid stated vector was introduced.')
            end

            obj.ElementType = myElementType;
        end

        % Normalize elements 
        function [obj] = Normalize(obj, direction, myLc)
            if (obj.Normalized && direction)
                obj.Normalized = false;
                obj = obj.Normalize(false, obj.Lc);
            end

            if (direction)
                obj.Lc = myLc;
                obj.Tc = sqrt(myLc^3/obj.mu);
                obj.mu = 1; 

                obj.Re = obj.Re/obj.Lc;

                switch (obj.ElementType)
                    case 'ECI'
                        obj.ElementSet = obj.ElementSet./[obj.Lc*ones(1,3) obj.Lc/obj.Tc*ones(1,3)];
                        obj.StateEvolution(:,2:end) = obj.StateEvolution(:,2:end)./[obj.Lc*ones(1,3) obj.Lc/obj.Tc*ones(1,3)];
                    case 'POL'
                        obj.ElementSet = obj.ElementSet([1 3 4 6])./[obj.Lc*ones(1,2) obj.Lc/obj.Tc*ones(1,2)];
                        obj.StateEvolution(:,2:end) = obj.StateEvolution(:,[1 3 4 6])./[obj.Lc*ones(1,2) obj.Lc/obj.Tc*ones(1,2)];                       
                    case 'COE'
                        obj.ElementSet([1 7]) = obj.ElementSet([1 7])/obj.Lc;
                        obj.StateEvolution(:,[2 obj.m+2]) = obj.StateEvolution(:,[2 obj.m+2])./obj.Lc;
                    case 'MOE'
                        obj.ElementSet(1) = obj.ElementSet(1)./obj.Lc;
                        obj.StateEvolution(:,2) = obj.StateEvolution(:,2)./obj.Lc;
                    otherwise
                        error('The current element set cannot be normalized.');
                end

                obj.StateEvolution(:,1) = obj.StateEvolution(:,1)/obj.Tc;

                obj.Normalized = true;
            elseif (obj.Normalized && ~direction)
                switch (obj.ElementType)
                    case 'ECI'
                        obj.ElementSet = obj.ElementSet.*[obj.Lc*ones(1,3) obj.Lc/obj.Tc*ones(1,3)];
                        obj.StateEvolution(:,2:obj.m+1) = obj.StateEvolution(:,2:obj.m+1).*[obj.Lc*ones(1,3) obj.Lc/obj.Tc*ones(1,3)];
                    case 'POL'
                        obj.ElementSet = obj.ElementSet([1 3 4 6]).*[obj.Lc*ones(1,2) obj.Lc/obj.Tc*ones(1,2)];
                        obj.StateEvolution(:,2:end) = obj.StateEvolution(:,[1 3 4 6]).*[obj.Lc*ones(1,2) obj.Lc/obj.Tc*ones(1,2)]; 
                    case 'COE'
                        obj.ElementSet([1 7]) = obj.ElementSet([1 7])*obj.Lc;
                        obj.StateEvolution(:,[2 obj.m+2]) = obj.StateEvolution(:,[2 obj.m+2])*obj.Lc;
                    case 'MOE'
                        obj.ElementSet(1) = obj.ElementSet(1).*obj.Lc;
                        obj.StateEvolution(:,2) = obj.StateEvolution(:,2)*obj.Lc;
                    otherwise
                        error('The current element set cannot be normalized.');
                end

                obj.StateEvolution(:,1) = obj.StateEvolution(:,1)*obj.Tc;

                obj.mu = (obj.Lc^3/obj.Tc^2);

                obj.Re = obj.Re*obj.Lc;

                obj.Normalized = false;
            end
        end
        
        % Plot the state evolution 
        function PlotTrajectory(obj, figureObj, InitialEpoch, FinalEpoch)
            % Sanity check 
            if (size(obj.StateEvolution,1) >= 2)
                if (FinalEpoch > obj.PropagatedEpoch)
                    FinalEpoch = obj.PropagatedEpoch;
                    warning('Propagation is not available spanning the selected epoch span.')
                end

                if (InitialEpoch > obj.PropagatedEpoch)
                    InitialEpoch = obj.Tspan(1);
                end
        
                % Plot the trajectory
                diff = abs(obj.Tspan-FinalEpoch); 
                [~, PosFinal] = sort(diff); 
                PosFinal = PosFinal(1);
    
                diff = abs(obj.Tspan-InitialEpoch); 
                [~, PosInit] = sort(diff);  
                PosInit = PosInit(1);
    
                % Change to Cartesian coordinates
                switch (obj.ElementType)
                    case 'ECI'
                    otherwise
                    obj = obj.State2Cartesian();
                end
    
                % Plot the Cartesian trajectory
                figureObj;
                view(3)
                plot3(obj.StateEvolution(PosInit:PosFinal,2), obj.StateEvolution(PosInit:PosFinal,3), obj.StateEvolution(PosInit:PosFinal,4));
                grid on; 
                xlabel('$x$');
                ylabel('$y$');
                zlabel('$z$');
            else
                warning('No propagation is available');
            end
        end
    end

    % Set graphics 
    methods (Static)
        function set_graphics()
            % Set graphical properties
            set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
            set(groot, 'defaultAxesFontSize', 11); 
            set(groot, 'defaultAxesGridAlpha', 0.3); 
            set(groot, 'defaultAxesLineWidth', 0.75);
            set(groot, 'defaultAxesXMinorTick', 'on');
            set(groot, 'defaultAxesYMinorTick', 'on');
            set(groot, 'defaultFigureRenderer', 'painters');
            set(groot, 'defaultLegendBox', 'off');
            set(groot, 'defaultLegendInterpreter', 'latex');
            set(groot, 'defaultLegendLocation', 'best');
            set(groot, 'defaultLineLineWidth', 1); 
            set(groot, 'defaultLineMarkerSize', 3);
            set(groot, 'defaultTextInterpreter','latex');
        end        
    end

    % Private methods 
    methods (Access = private)
        % Add the propagator model 
        function [f] = AddModel(obj, myModel)
            switch (myModel)
                case 'Keplerian'
                    f = @(tspan)obj.KeplerDynamics(tspan);
                case 'Osculating J2'
                    f = @(tspan)obj.APSDynamics(tspan, 0);
                case 'Mean J2'
                    f = @(tspan)obj.APSDynamics(tspan, 1);
                 case 'SGP4'
                    f = @(tspan)obj.APSDynamics(tspan, 2);
                case 'High-precision'
                    f = @(tspan)obj.NumericalProp(tspan);
                otherwise
                    error('The selected propagator model is not currently implemented.');
            end
        end

        % Change the state evolution to the Cartesian format
        function [obj] = State2Cartesian(obj)
            switch (obj.ElementType)
                case 'COE'
                    obj.ElementSet = obj.ECI2COE(obj.ElementSet, false);
                    aux = obj.ECI2COE(obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
                case 'POL'
                    obj.ElementSet = obj.ECI2POL(obj.ElementSet, false);
                    aux = obj.ECI2POL(obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];    
                case 'MOE'
                    obj.ElementSet = obj.ECI2MOE(obj.ElementSet, false);
                    obj.StateEvolution(:,2:end) = obj.ECI2MOE(obj.StateEvolution(:,2:end), false);
                case 'KS'
                    obj.ElementSet = obj.ECI2KS(obj.ElementSet, false);
                    obj.StateEvolution = [obj.StateEvolution(:,1) obj.ECI2KS(obj.StateEvolution(:,2:end), false)];
                otherwise
                    error('Requested transformation is not currently supported.');
            end
        end

        % Change the state evolution to the COE format
        function [obj] = State2COE(obj)
           switch (obj.ElementType)
                case 'ECI'
                    obj.ElementSet = obj.ECI2COE(obj.ElementSet, true);
                    aux = obj.ECI2COE(obj.StateEvolution(:,2:end), true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
                case 'POL'
                    obj.ElementSet = obj.ECI2POL(obj.ElementSet, false);
                    obj.ElementSet = obj.ECI2COE(obj.ElementSet, true);
                    aux = obj.ECI2POL(obj.StateEvolution(:,2:end), false);
                    aux = obj.ECI2COE(aux, true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
                case 'MOE'
                    obj.ElementSet = obj.COE2MOE(obj.ElementSet, false);
                    aux = obj.COE2MOE(obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
                case 'KS'
                    aux = obj.ECI2KS(obj.ElementSet, false);
                    obj.ElementSet = obj.ECI2COE(aux, true);
                    aux = obj.ECI2KS(obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution = [obj.StateEvolution(:,1) obj.ECI2COE(aux, true)];
                otherwise
                    error('Requested transformation is not currently supported.');
           end
        end

        % Change the state evolution to the MOE format
        function [obj] = State2MOE(obj)
           switch (obj.ElementType)
                case 'ECI'
                    obj.ElementSet = obj.ECI2MOE(obj.ElementSet, true);
                    aux = obj.ECI2MOE(obj.StateEvolution(:,2:end), true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
                case 'POL'
                    obj.ElementSet = obj.ECI2POL(obj.ElementSet, false);
                    obj.ElementSet = obj.ECI2MOE(obj.ElementSet, true);
                    aux = obj.ECI2POL(obj.StateEvolution(:,2:end), false);
                    aux = obj.ECI2MOE(aux, true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
                case 'COE'
                    obj.ElementSet = obj.COE2MOE(obj.ElementSet, true);
                    aux = obj.COE2MOE(obj.StateEvolution(:,2:end), true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
                case 'KS'
                    obj.ElementSet = obj.ECI2KS(obj.ElementSet, false);
                    obj.ElementSet = obj.ECI2COE(obj.ElementSet, true);
                    obj.ElementSet = obj.COE2MOE(obj.ElementSet, true);
                    aux = obj.ECI2KS(obj.StateEvolution(:,2:end), false);
                    aux = obj.ECI2COE(aux, true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) obj.COE2MOE(aux, true)];
                otherwise
                    error('Requested transformation is not currently supported.');
            end
        end

        % Change the state evolution to the POL format
        function [obj] = State2POL(obj)
           switch (obj.ElementType)
                case 'ECI'
                    obj.ElementSet = obj.ECI2POL(obj.ElementSet, true);
                    aux = obj.ECI2POL(obj.StateEvolution(:,2:end), true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
                case 'COE'
                    obj.ElementSet = obj.ECI2COE(obj.ElementSet, false);
                    obj.ElementSet = obj.ECI2POL(obj.ElementSet, true);
                    aux = obj.ECI2COE(obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution = [obj.StateEvolution(:,1) obj.ECI2POL(aux, true)];
                case 'MOE'
                    obj.ElementSet = obj.ECI2MOE(obj.ElementSet, false);
                    obj.ElementSet = obj.ECI2POL(obj.ElementSet, true);
                    aux = obj.ECI2MOE(obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution = [obj.StateEvolution(:,1) obj.ECI2POL(aux, true)];
                case 'KS'
                    obj.ElementSet = obj.ECI2KS(obj.ElementSet, false);
                    obj.ElementSet = obj.ECI2POL(obj.ElementSet, true);
                    aux = obj.ECI2KS(obj.StateEvolution(:,2:end), false);
                    aux = obj.ECI2POL(aux, true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
                otherwise
                    error('Requested transformation is not currently supported.');
            end
        end

        % Change the state evolution to the KS format
        function [obj] = State2KS(obj)
           lastwarn('The required transformation is not bijective.');
           lastwarn('Only for unperturbed circular orbits.');

           switch (obj.ElementType)
                case 'ECI'
                    obj.ElementSet = obj.ECI2KS(obj.ElementSet, true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) obj.ECI2KS(obj.StateEvolution(:,2:end), true)];
                case 'POL'
                    obj.ElementSet = obj.ECI2POL(obj.ElementSet, false);
                    obj.ElementSet = obj.ECI2KS(obj.ElementSet, true);
                    aux = obj.ECI2POL(obj.StateEvolution(:,2:end), false);
                    aux = obj.ECI2KS(aux, true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
                case 'COE'
                    obj.ElementSet = obj.ECI2COE(obj.ElementSet, false);
                    obj.ElementSet = obj.ECI2KS(obj.ElementSet, true);
                    aux = obj.ECI2COE(obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution = [obj.StateEvolution(:,1) obj.ECI2KS(aux, true)];
                case 'MOE'
                    obj.ElementSet = obj.ECI2MOE(obj.ElementSet, false);
                    obj.ElementSet = obj.ECI2KS(obj.ElementSet, true);
                    aux = obj.ECI2MOE(obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution = [obj.StateEvolution(:,1) obj.ECI2KS(aux, true)];
                otherwise
                    error('Requested transformation is not currently supported.');
            end
        end
    end
end




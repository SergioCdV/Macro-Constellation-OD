%% Constellation macro-orbit determination 
% Date: 01/30/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of an Orbit object

classdef Orbit
    % Fundamental definition 
    properties 
        mu                      % Central body gravitational parameter
    
        ElementType             % Type of orbital elements in use
        ElementSet              % Parametrizing orbital element set
        Sigma                   % Covariance matrix of the estimated element set

        InitialEpoch            % Initial epoch 
        CurrentEpoch            % Current propagation epoch 
        PropagatedEpoch;        % End time propagation epoch
        FinalEpoch              % Final orbit epoch
        TimeStep                % Integration time step

        StateEvolution;         % Evolution of the orbital element set
        Dynamics;               % Dynamical model to be used
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
                case 'Cartesian'
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
                    case 'Cartesian'
                        obj.ElementSet = myElementSet(1,1:obj.m);
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
            elseif (myCurrentEpoch < obj.FinalEpoch)
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

        % Propagate to current epoch 
        function [obj] = Propagate(obj)
            % Check time span limits 
            if (obj.PropagatedEpoch < obj.CurrentEpoch)
                aux_tspan = obj.PropagatedEpoch:obj.TimeStep:obj.CurrentEpoch;

                if (obj.Normalized)
                    aux_tspan = aux_tspan / obj.Tc;
                end
    
                % Perform the propagation 
                AuxEvolution = obj.Dynamics(aux_tspan); 

                % Update states
                obj.StateEvolution = [obj.StateEvolution; [obj.StateEvolution(end,1)+aux_tspan(2:end).' AuxEvolution(2:end,:)]];
                obj.ElementSet = AuxEvolution(end,:);

                obj.Tspan = aux_tspan;
    
                % Update the last propagated epoch 
                obj.PropagatedEpoch = obj.CurrentEpoch;
            end
        end

        % Change the state evolution
        function [obj] = ChangeStateFormat(obj, myElementType)
            switch (myElementType)
                case obj.ElementType
                case 'Cartesian'
                    obj = State2Cartesian(obj);
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

                switch (obj.ElementType)
                    case 'Cartesian'
                        obj.ElementSet = obj.ElementSet./[obj.Lc*ones(1,3) obj.Lc/obj.T*ones(1,3)];
                        obj.StateEvolution = obj.StateEvolution./[obj.Lc*ones(1,3) obj.Lc/obj.T*ones(1,3)];
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
                obj.Tspan = obj.Tspan/obj.Tc;

                obj.Normalized = true;
            else
                switch (obj.ElementType)
                    case 'Cartesian'
                        obj.ElementSet = obj.ElementSet.*[obj.Lc*ones(1,3) obj.Lc/obj.T*ones(1,3)];
                        obj.StateEvolution(:,2:obj.m+2) = obj.StateEvolution(:,2:obj.m+2)./[obj.Lc*ones(1,3) obj.Lc/obj.T*ones(1,3)];
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
                obj.Tspan = obj.Tspan*obj.Tc;

                obj.mu = (obj.Lc^3/obj.Tc^2);

                obj.Normalized = false;
            end
        end

        % Set graphics 
        function set_graphics(obj)
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

        % Plot the state evolution 
        function PlotTrajectory(obj, figureObj, InitialEpoch, FinalEpoch)
            % Sanity check 
            if (~isempty(obj.Tspan))
                if (FinalEpoch > obj.PropagatedEpoch)
                    FinalEpoch = obj.PropagatedEpoch;
                end
    
                if (FinalEpoch > obj.Tspan(1))
                    FinalEpoch = obj.Tspan(1);
                end
    
                % Plot the trajectory
                diff = abs(obj.Tspan-FinalEpoch); 
                [~, PosFinal] = sort(diff); 
                PosFinal = PosFinal(1);
    
                diff = abs(obj.Tspan-InitialEpoch); 
                [~, PosInit] = sort(diff);  
                PosInit = PosInit(1);
    
                % Change to Cartesian coordinates
                obj.Trajectory = obj.State2Cartesian(obj);
    
                % Plot the Cartesian trajectory
                figureObj;
                obj.set_graphics();
                plot3(obj.Trajectory(PosInit:PosFinal,1), obj.Trajectory(PosInit:PosFinal,2), obj.Trajectory(PosInit:PosFinal,3));
                grid on; 
                xlabel('$x$');
                ylabel('$y$');
                zlabel('$z$');
            else
                warning('No propagation is available');
            end
        end
    end

    % Private methods 
    methods (Access = private)
        % Add the propagator model 
        function [f] = AddModel(obj, myModel)
            switch (myModel)
                case 'Keplerian'
                    f = @(tspan)obj.KeplerDynamics(tspan);
                case 'J2'
                    f = obj.AddModel(obj, myModel);
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
                case 'Cartesian'
                    obj.ElementSet = obj.ECI2COE(obj.ElementSet, true);
                    aux = obj.ECI2COE(obj.StateEvolution(:,2:end), true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
                case 'MOE'
                    obj.ElementSet = obj.COE2MOE(obj.ElementSet, false);
                    aux = obj.COE2MOE(obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
                case 'KS'
                    aux = obj.ECI2KS(obj.ElementSet, false);
                    obj.ElementSet = obj.ECI2COE(aux, false);
                    aux = obj.ECI2KS(obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution = [obj.StateEvolution(:,1) obj.ECI2COE(aux, false)];
                otherwise
                    error('Requested transformation is not currently supported.');
           end
        end

        % Change the state evolution to the MOE format
        function [obj] = State2MOE(obj)
           switch (obj.ElementType)
                case 'Cartesian'
                    obj.ElementSet = obj.ECI2COE(obj.ElementSet, true);
                    obj.ElementSet = obj.COE2MOE(obj.ElementSet, true);
                    aux = obj.ECI2COE(obj.StateEvolution(:,2:end), true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) obj.COE2MOE(aux, true)];
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

        % Change the state evolution to the KS format
        function [obj] = State2KS(obj)
           lastwarn('The required transformation is not bijective.');

           switch (obj.ElementType)
                case 'Cartesian'
                    obj.ElementSet = obj.ECI2KS(obj.ElementSet, true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) obj.ECI2KS(obj.StateEvolution(:,2:end), true)];
                case 'COE'
                    obj.ElementSet = obj.ECI2COE(obj.ElementSet, false);
                    obj.ElementSet = obj.ECI2KS(obj.ElementSet, true);
                    aux = obj.ECI2COE(obj.StateEvolution(:,2:end), false);
                    aux = obj.ECI2KS(aux, true);
                    obj.StateEvolution = [obj.StateEvolution(:,1) aux];
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




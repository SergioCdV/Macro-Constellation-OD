%% Constellation macro-orbit determination 
% Date: 01/30/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of an Orbit object

classdef Orbit
    % Fundamental definition 
    properties 
        Mu                      % Central body gravitational parameter
    
        ElementType             % Type of orbital elements in use
        ElementSet              % Parametrizing orbital element set
        Sigma                   % Covariance matrix of the estimated element set

        InitialEpoch            % Initial epoch 
        CurrentEpoch            % Current propagation epoch 
        PropagatedEpoch;        % End time propagation epoch
        FinalEpoch              % Final orbit epoch
        TimeStep                % Integration time step

        StateEvolution;         % Evolution of the orbital element set
    end

    % Private properties 
    properties (Hidden = true)
        m = 6;                  % State space dimension

        Tspan                   % Integration time span

        IntegrationOptions;     % Integration tolerance

        Dynamics;               % Dynamical model to be used
        Trajectory;             % Cartesian trajectory

        Lc;                     % Characteristic length 
        Tc;                     % Characteristic time 
        Normalized = false;     % Flag to account for normalized coordinates
    end

    % Public methods 
    methods
        % General constructor
        function [obj] = Orbit(myElementType, myElementSet, myInitialEpoch)
            % User specified values
            obj.ElementSet = myElementSet;
            obj.InitialEpoch = myInitialEpoch;
            obj.CurrentEpoch = obj.InitialEpoch;

            % Dafault values 
            obj.IntegrationOptions = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      % Integration tolerances

            % Sanity checks 
            if (length(myElementSet) < obj.m)
                error('The state vector dimension is not valid.');
            end

            if (size(myElementSet,1) <= 1)
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

            obj.StateEvolution(1,:) = [myInitialEpoch, myElementSet];
        end

        % Add final epoch 
        function [obj] = AddFinalEpoch(obj, myFinalEpoch)
            if (myFinalEpoch < obj.InitialEpoch)
                error('Specified final epoch is not greater than the initial epoch');
            else
                obj.FinalEpoch = myFinalEpoch;
            end
        end

        % Specify current epoch 
        function [obj] = AddCurrentEpoch(obj, myCurrentEpoch)
            if (myCurrentEpoch < obj.CurrentEpoch)
                error('New current epoch is not greater than the actual current epoch');
            elseif (myCurrentEpoch < obj.FinalEpoch)
                obj.CurrentEpoch = myCurrentEpoch;
            elseif (myCurrentEpoch > obj.FinalEpoch)
                obj.FinalEpoch = myCurrentEpoch;
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
            switch (myModel)
                case 'Keplerian'
                    obj = obj.AddModel(obj, myModel);
                case 'J2'
                    obj = obj.AddModel(obj, myModel);
                otherwise
                    error('The selected propagator model is not currently implemented.');
            end
        end

        % Add integration tolerances
        function [obj] = AddIntegrationTolerances(myIntegrationOptions)
            obj.IntegrationOptions = myIntegrationOptions;
        end

        % Propagate to current epoch 
        function [obj] = Propagate(obj)
            % Check time span limits 
            if (obj.Tspan(end) < obj.CurrentEpoch)
                obj.Tspan = obj.InitialEpoch:obj.TimeStep:obj.CurrentEpoch;
            end

            % Perform the propagation 
            obj = obj.Dynamics(obj);  

            % Update the last propagated epoch 
            obj.PropagatedEpoch = obj.Tspan(end);
        end
        
        % Plot the state evolution 
        function PlotTrajectory(obj, InitialEpoch, FinalEpoch)
            % Sanity check 
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
            figure
            obj.set_graphics();
            plot3(obj.Trajectory(PosInit:PosFinal,1), obj.Trajectory(PosInit:PosFinal,2), obj.Trajectory(PosInit:PosFinal,3));
            grid on; 
            xlabel('$x$');
            ylabel('$y$');
            zlabel('$z$');
        end

        % Change the state evolution
        function [obj] = ChangeStateFormat(obj, myElementType)
            if (myElementType ~= obj.ElementType)
                switch (myElementType)
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
        end

        % Normalize elements 
        function [obj] = Normalize(obj, direction, myLc)
            if (obj.Normalized)
                obj = obj.Normalize(obj, obj.Lc, false);
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

                obj.Normalized = false;
            end
        end
    end

    % Private methods 
    methods (Access = private)
        % Add the propagator model 
        function [obj] = AddModel(obj, myModel)
            switch (myModel)
            end
        end

        % Change the state evolution to the Cartesian format
        function [obj] = State2Cartesian(obj)
            switch (obj.ElementType)
                case 'COE'
                    obj.ElementSet = ECI2COE(obj.mu, obj.ElementSet, false);
                    obj.StateEvolution(:,2:end) = ECI2COE(obj.mu, obj.StateEvolution(:,2:end), false);
                case 'MOE'
                    obj.ElementSet = ECI2MOE(obj.mu, obj.ElementSet, false);
                    obj.StateEvolution(:,2:end) = ECI2MOE(obj.mu, obj.StateEvolution(:,2:end), false);
                case 'KS'
                    obj.ElementSet = ECI2KS(obj.mu, obj.ElementSet, false);
                    obj.StateEvolution(:,2:end) = ECI2KS(obj.mu, obj.StateEvolution(:,2:end), false);
                otherwise
                    error('Requested transformation is not currently supported.');
            end
        end

        % Change the state evolution to the COE format
        function [obj] = State2COE(obj)
           switch (obj.ElementType)
                case 'Cartesian'
                    obj.ElementSet = ECI2COE(obj.mu, obj.ElementSet, true);
                    obj.StateEvolution(:,2:end) = ECI2COE(obj.mu, obj.StateEvolution(:,2:end), true);
                case 'MOE'
                    obj.ElementSet = COE2MOE(obj.mu, obj.ElementSet, false);
                    obj.StateEvolution(:,2:end) = COE2MOE(obj.mu, obj.StateEvolution(:,2:end), false);
                case 'KS'
                    aux = ECI2KS(obj.ElementSet, false);
                    obj.ElementSet = ECI2COE(obj.mu, aux, false);
                    aux = ECI2KS(obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution(:,2:end) = ECI2COE(obj.mu, aux, false);
                otherwise
                    error('Requested transformation is not currently supported.');
            end
        end

        % Change the state evolution to the MOE format
        function [obj] = State2MOE(obj)
           switch (obj.ElementType)
                case 'Cartesian'
                    obj.ElementSet = ECI2COE(obj.mu, obj.ElementSet, true);
                    obj.ElementSet = COE2MOE(obj.mu, obj.ElementSet, true);
                    obj.StateEvolution(:,2:end) = ECI2COE(obj.mu, obj.StateEvolution(:,2:end), true);
                    obj.StateEvolution(:,2:end) = COE2MOE(obj.mu, obj.StateEvolution(:,2:end), true);
                case 'COE'
                    obj.ElementSet = COE2MOE(obj.mu, obj.ElementSet, true);
                    obj.StateEvolution(:,2:end) = COE2MOE(obj.mu, obj.StateEvolution(:,2:end), true);
                case 'KS'
                    obj.ElementSet = ECI2KS(obj.mu, obj.ElementSet, false);
                    obj.ElementSet = ECI2COE(obj.mu, obj.ElementSet, true);
                    obj.ElementSet = COE2MOE(obj.mu, obj.ElementSet, true);
                    obj.StateEvolution(:,2:end) = ECI2KS(obj.mu, obj.StateEvolution(:,2:end), true);
                    obj.StateEvolution(:,2:end) = ECI2COE(obj.mu, obj.StateEvolution(:,2:end), true);
                    obj.StateEvolution(:,2:end) = COE2MOE(obj.mu, obj.StateEvolution(:,2:end), true);
                otherwise
                    error('Requested transformation is not currently supported.');
            end
        end

        % Change the state evolution to the KS format
        function [obj] = State2KS(obj)
           lastwarn('The required transformation is not bijective.');

           switch (obj.ElementType)
                case 'Cartesian'
                    obj.ElementSet = ECI2KS(obj.mu, obj.ElementSet, true);
                    obj.StateEvolution(:,2:end) = ECI2KS(obj.mu, obj.StateEvolution(:,2:end), true);
                case 'COE'
                    obj.ElementSet = ECI2COE(obj.mu, obj.ElementSet, false);
                    obj.ElementSet = ECI2KS(obj.mu, obj.ElementSet, true);
                    obj.StateEvolution(:,2:end) = ECI2COE(obj.mu, obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution(:,2:end) = ECI2KS(obj.mu, obj.StateEvolution(:,2:end), true);
                case 'MOE'
                    obj.ElementSet = ECI2MOE(obj.mu, obj.ElementSet, false);
                    obj.ElementSet = ECI2KS(obj.mu, obj.ElementSet, true);
                    obj.StateEvolution(:,2:end) = ECI2MOE(obj.mu, obj.StateEvolution(:,2:end), false);
                    obj.StateEvolution(:,2:end) = ECI2KS(obj.mu, obj.StateEvolution(:,2:end), true);
                otherwise
                    error('Requested transformation is not currently supported.');
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
    end
    
end




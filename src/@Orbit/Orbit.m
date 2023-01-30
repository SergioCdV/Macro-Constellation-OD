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

        StateEvolution = {};    % Evolution of the orbital element set
    end

    % Private properties 
    properties (Hidden = true)
        m = 6;                  % State space dimension

        Tspan                   % Integration time span

        IntegrationOptions;     % Integration tolerance

        Dynamics;               % Dynamical model to be used
        Trajectory;             % Cartesian trajectory
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
            if (length(myElementSet) ~= obj.m)
                error('The state vector dimension is not valid.');
            end

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

            obj.StateEvolution{1} = {myElementType, myInitialEpoch, myElementSet};
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
                obj.PropagatedEpoch = obj.Tspan(end);
            end

            % Perform the propagation 
            obj = obj.Dynamics(obj);  
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

        % Change the state evolution to the COE format

        % Change the state evolution to the MOE format

        % Change the state evolution to the KS format

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




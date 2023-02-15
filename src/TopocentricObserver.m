%% Constellation macro-orbit determination 
% Date: 15/02/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of an TopocentricObserver object

classdef TopocentricObserver
    properties
        Type = 'RADEC';         % Class of observer

        State = [];             % State of the observer
        Measurements = {};      % Measurements acquired by the observer

        InstrumentParams;       % Instrument parameters
        Sigma;                  % Covariance of the instrument

        InitialEpoch;           % Initial epoch
        PropagatedEpoch;        % Propgated epoch
        CurrentEpoch;           % Current epoch

        PD;                     % Probability of detection
        FOV = deg2rad(90);      % Field of view of the sensor
    end

    % Private properties
    properties (Access = private)
        StateDim = 2;           % Dimension of the state vector
        MeasDim;                % Dimension of the measurement vector
    end

    methods 
        % Constructor 
        function [obj] = TopocentricObserver(mySensor)
            switch (mySensor)
                case 'RADEC'
                    obj.Type = mySensor;
                    obj.MeasDim = 4; 
                case 'RADAR'
                    obj.Type = mySensor; 
                    obj.MeasDim = 2;
                otherwise 
                    error('No valid topocentric sensor has been selected.');
            end
        end

        % Add field of view 
        function [obj] = AddFOV(obj, myFOV)
            if (myFOV <= pi)
                obj.FOV = myFOV;
            else
                error('FOV cannot be greater than 180 deg.');
            end
        end

        % Add probability of detection 
        function [obj] = probability_detection(obj, myPD)
            obj.PD = myPD;
        end

        % Add instrument covariance 
        function [obj] = AddCovariance(obj, mySigma)
            try chol(mySigma);
                if (size(mySigma,1) == 1)
                    obj.Sigma = mySigma * eye(obj.StateDim);
                else
                    obj.Sigma = mySigma;
                end
            catch
                error('Input covariance matrix is not positive definite.'); 
            end
        end

        % Add an initial state 
        function [obj] = AddInitialState(obj, myInitialEpoch, myInitialState)
            if (myInitialEpoch >= 0)
                obj.InitialEpoch = myInitialEpoch; 
            end

            if (size(myInitialState,2) == obj.StateDim)
                obj.State = myInitialState;
            else
                error('No valid observer state has been addded.');
            end
        end

        % Observations
        function [timestamp, meas, StateEvolution] = Observe(obj, Orbit, Epoch)
            % Propagate the orbit 
            AuxOrbit = Orbit.SetCurrentEpoch(Epoch).Propagate().ChangeStateFormat('Cartesian');

            % Propagate the observer and take the measurements
            AuxOrbit = AuxOrbit.Normalize(false, 1);
            Tspan = AuxOrbit.StateEvolution(:,1);
            [Tspan, StateEvolution] = obj.Dynamics(AuxOrbit.InitialEpoch, obj.State, Tspan); 

            if (~isempty(StateEvolution))
                % Restrict the orbit state evolution   
                AuxOrbitEvolution = AuxOrbit.StateEvolution(Epoch + Tspan / (24 * 3600) >= obj.InitialEpoch,2:end);
                [timestamp, meas] = obj.ObservationProcess(Tspan, AuxOrbitEvolution, StateEvolution);
    
                % Apply the probability of detection
                index = logical(randsrc(length(timestamp), 1, [0, 1; 1-obj.PD, obj.PD]));
                timestamp = timestamp(index); 
                meas = meas(index,:); 
                StateEvolution = StateEvolution(index,:);
            else
                timestamp = []; 
                meas = []; 
                StateEvolution = [];
            end
        end

        % Observation process 
        function [t, meas] = ObservationProcess(obj, Tspan, Orbit, StateEvolution)
            % Preallocation
            meas = []; 
            t = []; 

            % Observation
            for i = 1:length(Tspan)
                if (dot(Orbit(i,1:3)/norm(Orbit(i,1:3)), StateEvolution(i,1:3)/norm(StateEvolution(i,1:3))) > cos(obj.FOV/2))
                    meas = [meas; obj.TopocentricObservation(obj.Type, StateEvolution(i,1:3), StateEvolution(i,4:6), Orbit(i,1:3), Orbit(i,4:6))];
                    t = [t; Tspan(i)];
                end
            end
        end

        % Topocentric observation
        function [meas] = TopocentricObservation(obj, sensor, Robs, Vobs, rsat, vsat)
            slant = (rsat-Robs).';                 % Slant vector
            vslant = (vsat-Vobs).';                % Slant rate vector 

            radar = norm(slant);                   % Radar slant measurement
            rate = dot(vslant, slant)/radar;       % Range rate measurement

            % Branch the algorithm
            switch (sensor)
                case 'RADAR'
                    meas(1,1) = radar;             % Radar slant measurement
                    meas(1,2) = rate;              % Range rate measurement
                    
                case 'RADEC'
                    if (sqrt(slant(1,1)^2+slant(2,1)^2) > 0)
                        sin_alpha = slant(2,1)/sqrt(slant(1,1)^2+slant(2,1)^2);
                        cos_alpha = slant(1,1)/sqrt(slant(1,1)^2+slant(2,1)^2);
                    else
                        sin_alpha = vslant(2,1)/sqrt(slant(1,1)^2+slant(2,1)^2);
                        cos_alpha = vslant(1,1)/sqrt(slant(1,1)^2+slant(2,1)^2);
                    end
 
                    meas(1,1) = atan2(sin_alpha, cos_alpha);   % Topoocentric right ascension
                    meas(1,3) = asin(slant(3,1)/radar(1,1));   % Topocentric declination

                    % Angle velocities
                    meas(1,2) = (vslant(2,1)*slant(1,1)-vslant(1,1)*slant(2,1))/(slant(1,1)^2+slant(2,1)^2);
                    meas(1,4) = (vslant(3,1)-rate*sin(meas(1,3)))/sqrt(slant(1,1)^2+slant(2,1)^2);

                otherwise 
                    error('No valid sensor has been selected.');
            end
        end

        % Gaussian likelihood function
        function [q] = LikelihoodFunction(obj, Sigma, z, y)
            res = y-z;
            q = exp((0.5*res.'*P^(-1)*res))/sqrt(det(Sigma)*(2*pi)^(size(Sigma,1)));
        end
    end
    
    methods (Access = private)
        % Dynamics 
        function [epoch, StateEvolution] = Dynamics(obj, Epoch, State, Tspan)
        % Check if an initial state has been given 
            if (~isempty(obj.State))
                % Constants 
                Re = 6378e3;        % Surface altitude over the ECI frame 
    
                % ECEF frame state
                r_ecef = Re * [cos(State(1,1)) * cos(State(1,2)) cos(State(1,1)) * sin(State(1,2)) sin(State(1,1))];
       
                % Synchronization 
                index = Epoch + Tspan / (24 * 3600) >= obj.InitialEpoch;
                epoch = Epoch + Tspan(index) / (24 * 3600);

                r_ecef = repmat(r_ecef, length(epoch), 1);
                v_ecef = zeros(length(epoch),3); 
                StateEvolution = zeros(length(epoch),6);
                
                utc = datetime(epoch, 'Format','yyyy MM dd HH mm ss.SSS', 'ConvertFrom', 'juliandate');
                utc.TimeZone = 'Z';
                utc = [year(utc) month(utc) day(utc) hour(utc) minute(utc) second(utc)];

                for i = 1:length(epoch)
                    [StateEvolution(i,1:3), StateEvolution(i,4:6)] = ecef2eci(utc(i,:), r_ecef(i,:), v_ecef(i,:));
                end
            else
                error('No obsever state has been given.');
            end
        end
    end
end
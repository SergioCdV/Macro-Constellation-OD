%% Constellation macro-orbit determination 
% Date: 15/02/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of an RadarSensor object

classdef RadarSensor < Sensors.AbstractSensor
    properties
        FOV = deg2rad(90);      % Field of view of the sensor
    end

    methods 
        % Constructor 
        function [obj] = RadarSensor(myInitialEpoch, myInitialState, mySigma, myPD, myFOV)
            myMeasDim = 2; 
            myStateDim = 2;
            obj = obj@Sensors.AbstractSensor(myStateDim, myMeasDim, myInitialEpoch, myInitialState, mySigma, myPD);

            if (exist('myFOV', 'var'))
                if (myFOV <= pi)
                    obj.FOV = myFOV;
                else
                    error('FOV cannot be greater than 180 deg.');
                end
            end
        end

        % Configuration of the clutter model 
        function [obj] = ConfigClutter(myPC, myNC)
            if (myPC >= 0 && myPC <= 1)
                obj.PC = myPC;
            else
                error('No valid detection probability has been input.'); 
            end

            if (myNC > 0)
                obj.NC = myPC; 
            else
                error('The average number of false measurements is not valid.');
            end
        end

        % Observations
        function [timestamp, meas, StateEvolution] = Observe(obj, Orbit, Epoch)
            % Propagate the orbit 
            AuxOrbit = Orbit.SetCurrentEpoch(Epoch).Propagate().ChangeStateFormat('ECI');

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

                % Add noise
                for i = 1:size(meas,1)
                    meas(i,:) = mvnrnd(meas(i,:), obj.Sigma, 1);
                end

                % Add clutter
                index = logical(randsrc(size(meas,1), 1, [0, 1; 1-obj.PC, obj.PC]));
                clutter = [normrnd(500e3, 5e2, length(index), 1) normrnd(7e3, 1e2, length(index), 1)];
                clutterTime = timestamp(index,:);
                clutterState = StateEvolution(index,:);
                
                if (size(clutter,1) > round(1.5 * obj.NC))
                    index = randi([0 size(clutter,1)], obj.NC, 1); 
                    clutter = clutter(index,:);
                    clutterTime = clutterTime(index,:);
                    clutterState = clutterState(index,:);
                end

                % Final assembly
                meas = [meas clutter(index,:)];
                timestamp = [timestamp; clutterTime];
                StateEvolution = [StateEvolution; clutterState];
                [timestamp, index] = sort(timestamp);
                meas = meas(index,:);
                StateEvolution = StateEvolution(index,:);
            else
                timestamp = []; 
                meas = []; 
                StateEvolution = [];
            end
        end

        % Topocentric observation
        function [meas] = TopocentricObservation(obj, Robs, Vobs, rsat, vsat)
            slant = (rsat-Robs).';                 % Slant vector
            vslant = (vsat-Vobs).';                % Slant rate vector 
            radar = norm(slant);                   % Radar slant measurement
            rate = dot(vslant, slant)/radar;       % Range rate measurement
            meas(1,1) = radar;                     % Radar slant measurement
            meas(1,2) = rate;                      % Range rate measurement
        end

        % Gaussian likelihood function
        function [q] = LikelihoodFunction(obj, Sigma, z, y)
            res = y-z;
            q = exp((0.5*res.'*P^(-1)*res))/sqrt(det(Sigma)*(2*pi)^(size(Sigma,1)));
        end
    end
    
    methods (Access = private)
        % Observation process 
        function [t, meas] = ObservationProcess(obj, Tspan, Orbit, StateEvolution)
            % Preallocation
            meas = []; 
            t = []; 

            % Observation
            for i = 1:length(Tspan)
                if (dot(Orbit(i,1:3)/norm(Orbit(i,1:3)), StateEvolution(i,1:3)/norm(StateEvolution(i,1:3))) < cos(obj.FOV/2))
                    meas = [meas; obj.TopocentricObservation(StateEvolution(i,1:3), StateEvolution(i,4:6), Orbit(i,1:3), Orbit(i,4:6))];
                    t = [t; Tspan(i)];
                end
            end
        end

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
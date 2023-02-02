%% Constellation macro-orbit determination 
% Date: 02/01/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of a Constellation object

classdef Constellation 
    % Fundamental properties 
    properties 
        Owner = 'MyNewSpaceStartUp'; 

        InitialEpoch = Inf;       % Initial establishment epoch of the constellation
        FinalEpoch = -Inf;        % Final epoch of the constellation
        CurrentEpoch = 0;         % Current epoch of the constellation

        TimeStep = 1e2;           % Sampling time step
            
        N = [];                   % Number of spacecraft (variable in time)
        Np = [];                  % Number of planes (variable in time)
        n = [];                   % Number of spacecraft per each plane (variable in plane and time)

        OrbitSet = {};            % Orbit set of the constellation (per spacecraft)
    end

    % Public methods
    methods 
        % Constructor
        function [obj] = Constellation(myConstellation, myN, myNp, myn)
            if (myN > 0 && myNp > 0 && myn > 0)
                switch (myConstellation)
                    case 'Walker'
                        obj = obj.Walker(myn * myNp, myNp, myn);
                    otherwise
                        obj.N = myN;
                        obj.Np = myNp; 
                        obj.n = myn; 
                end
            else
                lastwarn('Constellation definition makes no physical sense.');
            end
        end

        % Add an orbit to the orbital set
        function [obj] = AddOrbit(obj, Orbit)
            obj.OrbitSet = [obj.OrbitSet; {size(obj.OrbitSet,1)+1, Orbit}];

            if (Orbit.InitialEpoch > 0 && Orbit.InitialEpoch < obj.InitialEpoch)
                obj.InitialEpoch = Orbit.InitialEpoch; 
            end

            if (Orbit.FinalEpoch > 0 && Orbit.FinalEpoch > obj.FinalEpoch)
                obj.FinalEpoch = Orbit.FinalEpoch; 
            end
        end
        
        % Change time step 
        function [obj] = ChangeTimeStep(myTimeStep)
            if (myTimeStep > 0)
                obj.TimeStep = myTimeStep;
            end
        end

        % Remove an orbit my index
        function [obj] = RemoveOrbit(obj, OrbitID)
            obj.OrbitSet = obj.OrbitSet([1:OrbitID-1 OrbitID+1:end],:);
        end

        % Compute the number of planes in time 
        function [N, Np] = NumberOfPlanes(obj)
            if (~isempty(obj.OrbitSet))
                Tspan = obj.EpochSpan();
                N = zeros(1,length(Tspan)); 
                PlaneSet = cell(1,length(Tspan));
 
                % Number of planes
                for i = 1:length(Tspan)
                    PlaneSet{i} = zeros(1,3);
                    for j = 1:size(obj.OrbitSet,1)
                        Aux = obj.OrbitSet{j,2}.ChangeStateFormat('COE');
                        PlaneIndex = ismember(PlaneSet{i}, Aux.ElementSet(3:5), 'rows');
                        if (sum(PlaneIndex) == 0)
                            N(i) = N(i)+1;
                            PlaneSet{i} = [PlaneSet{i}; Aux.ElementSet(3:5)];
                        end
                    end
                end

                % Number of spacecraft per plane
                Np = zeros(length(Tspan), max(N));
                for i = 1:length(Tspan)
                    planes = PlaneSet{i}(2:end,:);
                    for j = 1:size(planes,1)
                        for k = 1:size(obj.OrbitSet,1)
                            Aux = obj.OrbitSet{k,2}.ChangeStateFormat('COE');
                            if (all(planes(j,:) == Aux.ElementSet(3:5)))
                                Np(i,j) = Np(i,j)+1;
                            end
                        end
                    end
                end
            else
                error('Orbital set of the constellation is empty.');
            end
        end

        % Compute the number of spacecraft in time
        function [n] = NumberOfSpacecraft(obj)
            if (~isempty(obj.OrbitSet))
                Tspan = obj.InitialEpoch:obj.TimeStep:obj.FinalEpoch;
                n = zeros(1,length(Tspan)); 

                for i = 1:length(Tspan)
                    for j = 1:size(obj.OrbitSet,1)
                        if (Tspan(i) <= obj.OrbitSet{j,2}.FinalEpoch && obj.OrbitSet{j,2}.InitialEpoch <= Tspan(i))
                            n(i) = n(i)+1;
                        end
                    end
                end
            else
                error('Orbital set of the constellation is empty.');
            end
        end

        % Propagate to a given epoch the set of orbits 
        function [obj] = Propagate(obj, myCurrentEpoch)
            if (myCurrentEpoch > 0 && myCurrentEpoch > obj.CurrentEpoch)
                % Sanity checks
                if (myCurrentEpoch > obj.FinalEpoch)
                    myCurrentEpoch = obj.FinalEpoch;
                end

                if (myCurrentEpoch < obj.InitialEpoch)
                    myCurrentEpoch = obj.InitialEpoch;
                end

                % Orbit propagation 
                for i = 1:size(obj.OrbitSet,1)
                    AuxOrbit = obj.OrbitSet{i,2}.SetCurrentEpoch(myCurrentEpoch);
                    obj.OrbitSet{i,2} = AuxOrbit.Propagate();
                end

                % Update the propagation epoch 
                obj.CurrentEpoch = myCurrentEpoch;
            end
        end
    end

    % Private methods
    methods (Access = private)
        function [tspan] = EpochSpan(obj)
            tspan = obj.InitialEpoch:obj.TimeStep:obj.FinalEpoch;
        end
    end
end
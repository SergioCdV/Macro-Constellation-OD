%% Constellation macro-orbit determination 
% Date: 02/01/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of a Constellation object

classdef Constellation 
    % Fundamental properties 
    properties 
        Owner = 'MyNewSpaceStartUp'; 

        InitialEpoch;       % Initial establishment epoch of the constellation
        FinalEpoch;         % Final epoch of the constellation
        CurrentEpoch;       % Current epoch of the constellation

        TimeStep = 1e-3;    % Sampling time step
            
        N = [];             % Number of spacecraft (variable in time)
        Np = [];            % Number of planes (variable in time)
        n = [];             % Number of spacecraft per each plane (variable in plane and time)

        OrbitSet = {};      % Orbit set of the constellation (per spacecraft)
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
        function [obj] = AddOrbit(obj, EpochSpan, Orbit)
            obj.OrbitSet = [obj.OrbitSet; {size(obj.OrbitSet,1)+1, EpochSpan, Orbit}];

            if (EpochSpan(1) > 0 && EpochSpan(1) < obj.InitialEpoch)
                obj.InitialEpoch = EpochSpan(1); 
            end

            if (EpochSpan(end) > 0 && EpochSpan(end) > obj.FinalEpoch)
                obj.FinalEpoch = EpochSpan(end); 
            end
        end

        function [obj] = RemoveOrbit(obj, OrbitID)
            obj.OrbitSet = obj.OrbitSet([1:OrbitID-1 OrbitID+1:end]);
        end

        % Compute the number of planes in time 
        function [N] = NumberOfPlanes(obj)
            
        end

        % Compute the number of spacecraft in time
        function [n] = NumberOfSpacecraft(obj)
            for i = 1:size(obj.OrbitSet,1)
            end
        end

        % Propagate to a given epoch the set of orbits 
        function []
    end

    % Private methods
    methods (Access = private)
        function [tspan] = EpochSpan(obj)
            tspan = obj.InitialEpoch:obj.TimeStep:obj.FinalEpoch;
        end
    end
end
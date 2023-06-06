

classdef (Abstract) Planet
    properties 
        Name = 'Earth';
        mu = 3.986004418e14;                % Gravitational Parameter 
        Re = 6371.137e3;                    % Mean radius
        omega = (2*pi/(24*3600));           % Angular velocity
        rho = [1.225 3.725e-12];            % Atmospheric density references
        mum = 7.96e15;                      % Magnetic field dipole strength 
    end
end
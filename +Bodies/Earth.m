


classdef Earth < Bodies.Planet

    properties 
        J2;

        erad;
        prad; 
        erot;
    end

    methods 
        function [obj] = Earth()
            % Constants 
            obj.Name = 'Earth';
            obj.mu = 3.986004418e14;                % Gravitational Parameter of the Earth
            obj.Re = 6371.137e3;                    % Mean Earth Radius
            obj.omega = (2*pi/(24*3600));           % Earth angular velocity
            obj.rho(1) = 1.225;                     % Air density at sea level
            obj.rho(2) = 3.725e-12;                 % Air density at 600 km
            obj.mum = 7.96e15;                      % Earth magnetic field dipole strength    
            obj.J2 = 1.08263e-3;                    % Second zonal harmonic of the Earth

            % Mean spherical earth
            obj.erad    = 6371008.7714;    % equatorial radius (meters)
            obj.prad    = 6371008.7714;    % polar radius (meters)
            obj.erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
        end

        function plot(obj)
            npanels = 180;                  % Number of globe panels around the equator deg/panel = 360/npanels
            alpha   = 0.25;                 % globe transparency level, 1 = opaque, through 0 = invisible
            GMST0 = 4.89496121282306;       % Set up a rotatable globe at J2000.0

            % Set initial view
            view(3);

            [x, y, z] = ellipsoid(0, 0, 0, obj.erad, obj.erad, obj.prad, npanels);
            globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
            
            if false
                hgx = hgtransform;
                set(hgx,'Matrix', makehgtform('zrotate',GMST0));
                set(globe,'Parent',hgx);
            end
           
            cdata = imread(obj.texture);  % Load Earth image for texture map
            set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
        end
    end
end
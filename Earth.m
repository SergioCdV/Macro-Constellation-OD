classdef Earth

    properties
        texture = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';

        % Mean spherical earth
        erad    = 6371008.7714; % equatorial radius (meters)
        prad    = 6371008.7714; % polar radius (meters)
        erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
    end

    methods 
        function plot(obj)
            npanels = 180;                  % Number of globe panels around the equator deg/panel = 360/npanels
            alpha   = 0.25;                    % globe transparency level, 1 = opaque, through 0 = invisible
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
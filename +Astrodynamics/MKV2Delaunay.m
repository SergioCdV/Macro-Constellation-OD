function [s] = MKV2Delaunay(mu, x, direction)
    if (direction)
        % Transformation from Delaunay to COE elements 
        s = Astrodynamics.MKV2COE(mu, x, true);

        % Transformation from COE to MKV
        s = Astrodynamics.Delaunay2COE(mu, s, false);
    else
        % Transformation from Delaunay to COE elements 
        s = Astrodynamics.Delaunay2COE(mu, x, true);

        % Transformation from COE to MKV
        s = Astrodynamics.MKV2COE(mu, s, false);
    end
end
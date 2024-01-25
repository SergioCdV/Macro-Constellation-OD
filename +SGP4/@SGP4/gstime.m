%% gstime function
% This function computes the Greenwich Mean Sidereal Time of a given epoch in Julian format 
% Author: Sergio Cuevas
% Created: 2023-12-11

% Inputs: - jdut1, the Julian date of interest
% Output: - double gst, the associated GMST

function [gst] = gstime(jdut1)
    tut1 = (jdut1 - 2451545.0) / 36525.0;       % Julian century
    temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + (876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841;  % sec
    temp = SGP4.SGP4.fmod(temp * SGP4.SGP4.deg2rad / 240.0, SGP4.SGP4.twopi);          % 360/86400 = 1/240, to deg, to rad

    % ------------------------ check quadrants ---------------------
    if (temp < 0.0)
        temp = temp + SGP4.SGP4.twopi;
    end

    gst = temp;
end 
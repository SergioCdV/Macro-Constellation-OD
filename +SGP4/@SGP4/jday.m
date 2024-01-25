%% jday function
% This function computes the Julian date associated to a given Gregorian epoch 
% Author: Sergio Cuevas
% Created: 2023-12-11

% Inputs: - int year, mon, day, hr, minute, the year, month, day, hour and
%           minute associated to the input Gregorian epoch
%         - double sec, the remaining seconds and fractino of seconds of
%           the epoch
% Output: - array output, with the Julian day and fraction of the day
%           associated to the input epoch

function [output] = jday( year, mon, day, hr, minute, sec)

    jd = 367.0 * year - floor((7 * (year + floor((mon + 9) / 12.0))) * 0.25) + floor(275 * mon / 9.0) + day + 1721013.5;  % use - 678987.0 to go to mjd directly
    jdFrac = (sec + minute * 60.0 + hr * 3600.0) / 86400.0;

    % Check that the day and fractional day are correct
    if (abs(jdFrac) > 1.0)
        dtt = floor(jdFrac);
        jd = jd + dtt;
        jdFrac = jdFrac - dtt;
    end

    output = [jd jdFrac];
end 
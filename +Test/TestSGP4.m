%% TestSGP4 
% This script provides verification and validation tests for the
% performance of the implemented SGP4 model. The source code may be found
% in 

% In particular, Vallado's set of test TLE are propagated through this
% matlab implementation and the output compared against the standard
% results from the original code. 3sigma bounds are also provided for both
% the position and velocity error metrics. 

% Please note that the given tests only verify the implementation of the
% SGP4 in Matlab against its oginal C code, but by no means the error
% metrics correspond to the accuray of the propagation itself.

clear 
clc 
close

%% Input data
file_name = '+Test\SGP4-VER.TLE';     % Test data file name
output_name = '+Test\tcppver.out';    % Output verification data

fid = fopen(file_name);         % Open file
tlecnt = 1;                     % TLE counter
tles = cell(tlecnt);            % Preallocation of a TLE cell array structure
line1 = fgetl(fid);

while ( ischar(line1) )
    if (line1(1) == '1')
        line2 = fgetl(fid);
        tle = SGP4.TLE(line1, line2);    % Create the TLE object
        tles{tlecnt} = tle;              % Save the object
        tlecnt = tlecnt + 1;
    end
    line1 = fgetl(fid);
end

fclose(fid);                        % Close file

%% Verification and validation of the propagator
curobj = 0;     % Current object index
strind = 0;
prevrv = 0;
ind = 0;
cnt = 0;

% Preallocation of the error metrics
rerr(1) = 0;
verr(1) = 0;

% Loop over the output verification file to compute Cartesian elements for various TLEs and epochs
fid = fopen(output_name);
line1 = fgetl(fid);

while ( ischar(line1) )
    strind = strfind(line1, 'xx'); 

    if (strind > 0)
        % Get the object number from the output file
        curobj = str2double(strtrim(line1(1:strind-1)));

        % Find the TLE with the matching object number
        ind = 1;
        while (ind <= tlecnt && tles{ind}.objectNum ~= curobj)
            ind = ind+1;
        end
        tle = tles{ind};

    elseif (curobj > 0)
        % Get the RV for the specified number of minutes (propagation)
        verout = sscanf(line1, '%f');
        rv = tle.getRV(verout(1));
        
        % Original code reused RV so if there is an error we need to carry old values
        if(tle.sgp4Error > 0)
            rv = prevrv;
        end

        cnt = cnt + 1;
        rdist = norm(rv(:,1)-verout(2:4));
        vdist = norm(rv(:,2)-verout(5:7));

        rerr(cnt) = rdist;
        verr(cnt) = vdist;

        if (rdist > 1e-7 || vdist > 1e-9)
            fprintf("Propagation error detected: \n")
            sprintf('Error metrics: %d\t%d\t%.15f\t%.15f\n', curobj, verout(1), rdist, vdist)
        end

        prevrv = rv;
    end

    line1 = fgetl(fid);
end
fclose(fid);

%% Analysis
mu_rerr = mean(1e6 * rerr);        % Mean position error [mm]
mu_verr = mean(1e6 * verr);        % Mean velocity error [mm]

sigma_rerr = std(1e6 * rerr);      % Standard deviation of the position error [mm/s]
sigma_verr = std(1e6 * verr);      % Standard deviation of the velocity error [mm/s]

%% Display
fprintf('Mean position error: %E mm\n',   (mu_rerr));
fprintf('3sigma position error: %E mm\n', (sigma_rerr));
fprintf('Mean velocity error: %E mm/s\n', (mu_verr));
fprintf('3sigma velocity error: %E mm/s\n', (sigma_verr));

figure
subplot(2,1,1)
hold on
plot(1:cnt, 1e6 * rerr + 3 * sigma_rerr)
yline(mu_rerr, 'k--')
scatter(1:cnt, 1e6 * rerr, 'filled')
ylabel('Position error [mm]')
legend('3sigma', 'Mean', 'Error')
grid on;

subplot(2,1,2)
hold on
plot(1:cnt, 1e6 * verr + 3 * sigma_verr)
yline(mu_verr, 'k--')
scatter(1:cnt, 1e6 * verr, 'filled')
ylabel('Velocity error [mm/s]')
legend('3sigma', 'Mean', 'Error')
grid on;

xlabel('TLE case')



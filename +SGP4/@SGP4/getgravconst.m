%% getgravconst function
% This function computes the constants of the Earth's reference system to be used 
% Author: Sergio Cuevas
% Created: 2023-12-11

% Inputs: - selector whichconst, which defines the model to be used in the SGP4 propagator
% Output: - object rec, the spacecraft propagation object

function getgravconst(whichconst, rec)
    rec.whichconst = whichconst;

    % -- wgs-72 low precisionconstants -- %
    if (whichconst == SGP4.SGP4().wgs72old)
        rec.mu = 398600.79964;            % in km3 / s2
        rec.radiusearthkm = 6378.135;     % km
        rec.xke = 0.0743669161;           % reciprocal of tumin
        rec.tumin = 1.0 / rec.xke;
        rec.j2 = 0.001082616;
        rec.j3 = -0.00000253881;
        rec.j4 = -0.00000165597;
        rec.j3oj2 = rec.j3 / rec.j2;

        % ------------ wgs-72 constants ------------ %
    elseif (whichconst == SGP4.SGP4().wgs72)
        rec.mu = 398600.8;                % in km3 / s2
        rec.radiusearthkm = 6378.135;     % km
        rec.xke = 60.0 / sqrt(rec.radiusearthkm*rec.radiusearthkm*rec.radiusearthkm / rec.mu);
        rec.tumin = 1.0 / rec.xke;
        rec.j2 = 0.001082616;
        rec.j3 = -0.00000253881;
        rec.j4 = -0.00000165597;
        rec.j3oj2 = rec.j3 / rec.j2;

    else 
        % ------------ wgs-84 constants ------------ %
        rec.mu = 398600.5;                % in km3 / s2
        rec.radiusearthkm = 6378.137;     % km
        rec.xke = 60.0 / sqrt(rec.radiusearthkm*rec.radiusearthkm*rec.radiusearthkm / rec.mu);
        rec.tumin = 1.0 / rec.xke;
        rec.j2 = 0.00108262998905;
        rec.j3 = -0.00000253215306;
        rec.j4 = -0.00000161098761;
        rec.j3oj2 = rec.j3 / rec.j2;
    end
end
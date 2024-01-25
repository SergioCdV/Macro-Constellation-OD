% octave/matlab version of SGP4 from Vallado
%
%
%     ----------------------------------------------------------------
%
%                               sgp4unit.cpp
%
%    this file contains the sgp4 procedures for analytical propagation
%    of a satellite. the code was originally released in the 1980 and 1986
%    spacetrack papers. a detailed discussion of the theory and history
%    may be found in the 2006 aiaa paper by vallado, crawford, hujsak,
%    and kelso.
%
%                            companion code for
%               fundamentals of astrodynamics and applications
%                                    2013
%                              by david vallado
%
%     (w) 719-573-2600, email dvallado@agi.com, davallado@gmail.com
%
%    current :
%               7 dec 15  david vallado
%                           fix jd, jdfrac
%

classdef SGP4

    properties (Constant)
        twopi = 2.0 * pi;       % 2pi
        deg2rad = pi / 180.0;   % DEG to RAD conversion factor
        wgs72old = 1;           % WGS72 model selector
        wgs72 = 2;              % WGS72 model selector
        wgs84 = 3;              % WGS84 model selector
    end
    
    methods (Static)
        % Deep space initialization functions
        dpper(e3, ee2, peo, pgho, pho, pinco, plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, ...
              t, xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4, zmol, zmos,init, rec, opsmode);
        
        dscom(epoch, ep, argpp, tc, inclp, nodep, np, rec);
        
        dsinit (tc, xpidot, rec);
        
        dspace(tc, rec);
        
        initl(epoch, rec);
        
        % Propagation per se
        [ok] = sgp4init( opsmode, satrec );
        [ok] = sgp4_propagator(satrec, tsince, rvhand)
        
        % Compute the mod operation with doubles
        function [fm] = fmod(numer, denom)
            tquot = floor(numer/denom);
            fm = numer-tquot * denom;
        end
            
        % Compute the Julian day 
        [output] = jday( year, mon, day, hr, minute, sec); 
        
        % Compute GMST of a Julian epoch
        [gst] = gstime(jdut1);
        
        % Get the Earth's reference frame constants
        getgravconst(whichconst, rec);
    end
end

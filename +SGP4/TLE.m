%% TLE class
% This class implements the TLE of a given spacecraft 
% Author: aholinch
% Created: 2018-09-09

classdef TLE < handle
  
    properties
        rec = 0;
        
        line1 = '';         % First line
        line2 = '';         % Second line
        
        intlid = '';        % Line number
        objectNum = 0;      % Number of the object
        epoch = 0;          % Initial epoch
        ndot = 0;           % First derivative of the mean motion
        nddot = 0;          % Second derivative of the mean motion
        bstar = 0;          % Bstar coefficient
        elnum = 0;
        incDeg = 0;         % Orbital inclination in deg
        raanDeg = 0;        % RAAN in deg
        ecc = 0;            % Orbital eccentricity
        argpDeg = 0;        % AoP in deg
        maDeg = 0;          % Mean anomaly in deg   
        n = 0;              % Mean motion
        revnum = 0;         % Number of revolutions
        sgp4Error = 0;      % Error code
    end

    methods    
        % Class constructor
        function obj = TLE(l1, l2)
            parseLines(obj,l1,l2);
        end
        
        % Function to parse the given TLE in the required information
        function parseLines(this, l1, l2)
	        this.rec = SGP4.elsetrec();
	        this.line1 = l1;
	        this.line2 = l2;
        
            % Object number
            this.objectNum = str2double(this.line1(3:7));
        
            % Derivative of the mean motion (floating point assumed)
            this.intlid = strtrim(this.line1(10:18));
            tmp = this.line1(34:43);
            tmp(1) = '0';
            this.ndot = str2double(tmp);
        
            % Second derivative of the mean motion (floating point assumed)
            tm = 1;
            if strcmp('-', this.line1(45:45))
                tm = -1;
            end
            tmp = this.line1(44:50);
            tmp(1) ='0';
            tmp(2) ='.';
            this.nddot = str2double(tmp);
            ev = str2double(this.line1(51:52));
            this.nddot = this.nddot * (10^ev);
            this.nddot = tm*this.nddot;
        
            % Bstar coefficient (floating point assumed)
            tm = 1;
            if strcmp('-', this.line1(54:54))
                tm = -1;
            end
            tmp = this.line1(53:59);
            tmp(1) ='0';
            tmp(2) ='.';
            this.bstar = str2double(tmp);
            ev = str2double(this.line1(60:61));
            this.bstar = this.bstar * (10^ev);
            this.bstar = 0 * tm * this.bstar;
        
            % Epoch 
            this.elnum = str2double(strtrim(this.line1(65:68)));
            this.epoch = parseEpoch(this,this.line1(19:32));
        
            % Eccentricity (floating point assumed)
            tmp = this.line2(25:33);
            tmp(1) ='0';
            tmp(2) ='.';
            this.ecc =           str2double(strtrim(tmp));
        
            % Remaining orbital elements
            this.incDeg =        str2double(strtrim(this.line2(9:16)));
            this.raanDeg =       str2double(strtrim(this.line2(18:25)));
            this.argpDeg =       str2double(strtrim(this.line2(35:42)));
            this.maDeg =         str2double(strtrim(this.line2(44:51)));
            this.n =             str2double(strtrim(this.line2(53:63)));
            this.revnum =        str2double(strtrim(this.line2(64:68)));
            this.rec.classification =               this.line1(8);
        
            % Transform to the correct units
            setValsToRec(this)
        end
        
        % Function to transform the parsed line to the correct units
        function setValsToRec(this)
            deg2rad = pi / 180.0;           % DEG to RAD conversion
            xpdotp = 1440.0 / (2.0 * pi);   % 229.1831180523293 
        
            this.rec.elnum = this.elnum;
            this.rec.revnum = this.revnum;
            this.rec.satnum = this.objectNum;
            this.rec.bstar =  0 * this.bstar;
            this.rec.inclo =  this.incDeg  * deg2rad;
            this.rec.nodeo =  this.raanDeg * deg2rad;
            this.rec.argpo =  this.argpDeg * deg2rad;
            this.rec.mo =     this.maDeg   * deg2rad;
            this.rec.ecco =   this.ecc;
        
            this.rec.no_kozai = this.n/xpdotp;
            this.rec.ndot =     this.ndot / (xpdotp*1440.0);
            this.rec.nddot =    this.nddot / (xpdotp*1440.0*1440.0);
        
            % Initialize the remaining temporal variables
            SGP4.SGP4().sgp4init('a', this.rec);
        end
        
        % Parse the TLE epoch
        function [ep] = parseEpoch(this, str)
            % Parse the TLE year
            year = str2double(strtrim(str(1:2)));
            if (year > 56)
                year = year + 1900;
            else
                year = year + 2000;
            end 
        
            % Parse the remaining epoch in Gregorian format (fraction of day, floating point assumed)
            doy = str2double(strtrim(str(3:5)));
            str(5) = '0';
            dfrac = str2double(strtrim(str(5:14)));
        
            dfrac = 24.0 * dfrac;       % Fraction of day in hours
            hr = floor(dfrac);          % Hour
            dfrac = 60.0 * (dfrac-hr);  % Fraction of hour
            min = floor(dfrac);         % Minutes
            sec = 60.0 * (dfrac-min);   % Seconds
        
            mons = [31 28 31 30 31 30 31 31 30 31 30 31];
        
            % Check if the epoch was a leap year
            il = ((rem (year, 4) == 0 & rem (year, 100) ~= 0) | rem (year, 400) == 0);
            if(il)
                mons(2) = 29;
            end
        
            mon = 1;
            days = doy;
            while( mon < 12 && days > mons(mon))
                days = days - mons(mon);
                mon = mon + 1;
            end
        
            % Compute the Julian date of the epoch
            jda = SGP4.SGP4.jday(year, mon, days, hr, min, sec);
        
            this.rec.jdsatepoch = jda(1);     % Julian day 
            this.rec.jdsatepochF = jda(2);    % Fraction of the year
        
            % Error code
            ep = 0;                         
        end
        
        % Propagate the spacecraft to a given elapsed epoch
        function [rv] = getRV(this, elapsed_time)
            % Initialize the proapagation
	        this.rec.error = 0;   
            rvhand = SGP4.rvhandle();
        
            % Propagate the spacecraft
	        SGP4.SGP4.sgp4_propagator(this.rec, elapsed_time, rvhand);
	        this.sgp4Error = this.rec.error;
        
            r = rvhand.r;   % Propagated position
            v = rvhand.v;   % Propagated velocity
            rv = [r v];     % Final TEME Cartesian elements
        end
    end
end

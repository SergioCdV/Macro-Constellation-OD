function dspace(tc, rec)    
    xndt = 0;
    xnddt = 0;
    xldot = 0;
    
    fasx2 = 0.13130908;
    fasx4 = 2.8843198;
    fasx6 = 0.37448087;
    g22 = 5.7686396;
    g32 = 0.95240898;
    g44 = 1.8014998;
    g52 = 1.0508330;
    g54 = 4.4108898;
    rptim = 4.37526908801129966e-3; % this equates to 7.29211514668855e-5 rad/sec
    stepp = 720.0;
    stepn = -720.0;
    step2 = 259200.0;

    %/* ----------- calculate deep space resonance effects ----------- */
    rec.dndt = 0.0;
    theta = SGP4.SGP4.fmod(rec.gsto + tc * rptim, SGP4.SGP4.twopi);
    rec.em = rec.em + rec.dedt * rec.t;

    rec.inclm = rec.inclm + rec.didt * rec.t;
    rec.argpm = rec.argpm + rec.domdt * rec.t;
    rec.nodem = rec.nodem + rec.dnodt * rec.t;
    rec.mm = rec.mm + rec.dmdt * rec.t;

    %   sgp4fix for negative inclinations
    %   the following if statement should be commented out
    %  if (inclm < 0.0)
    % {
    %    inclm = -inclm;
    %    argpm = argpm - pi;
    %    nodem = nodem + pi;
    %  }

    %/* - update resonances : numerical (euler-maclaurin) integration - */
    %/* ------------------------- epoch restart ----------------------  */
    %   sgp4fix for propagator problems
    %   the following integration works for negative time steps and periods
    %   the specific changes are unknown because the original code was so convoluted

    % sgp4fix take out atime = 0.0 and fix for faster operation
    ft = 0.0;
    if (rec.irez ~= 0)
        % sgp4fix streamline check
        if ((rec.atime == 0.0) || (rec.t * rec.atime <= 0.0) || (abs(rec.t) < abs(rec.atime)))
            rec.atime = 0.0;
            rec.xni = rec.no_unkozai;
            rec.xli = rec.xlamo;
        end

        % sgp4fix move check outside loop
        if (rec.t > 0.0)
            delt = stepp;
        else
            delt = stepn;
        end

        iretn = 381; % added for do loop
        while (iretn == 381)
            %/* ------------------- dot terms calculated ------------- */
            %/* ----------- near - synchronous resonance terms ------- */
            if (rec.irez ~= 2)
                xndt = rec.del1 * sin(rec.xli - fasx2) + rec.del2 * sin(2.0 * (rec.xli - fasx4)) + rec.del3 * sin(3.0 * (rec.xli - fasx6));
                xldot = rec.xni + rec.xfact;
                xnddt = rec.del1 * cos(rec.xli - fasx2) + 2.0 * rec.del2 * cos(2.0 * (rec.xli - fasx4)) + 3.0 * rec.del3 * cos(3.0 * (rec.xli - fasx6));
                xnddt = xnddt * xldot;
            else
                %/* --------- near - half-day resonance terms -------- */
                xomi = rec.argpo + rec.argpdot * rec.atime;
                x2omi = xomi + xomi;
                x2li = rec.xli + rec.xli;
                xndt = rec.d2201 * sin(x2omi + rec.xli - g22) + rec.d2211 * sin(rec.xli - g22) + ...
                        rec.d3210 * sin(xomi + rec.xli - g32) + rec.d3222 * sin(-xomi + rec.xli - g32) + ...
                        rec.d4410 * sin(x2omi + x2li - g44) + rec.d4422 * sin(x2li - g44) + ...
                        rec.d5220 * sin(xomi + rec.xli - g52) + rec.d5232 * sin(-xomi + rec.xli - g52) + ...
                        rec.d5421 * sin(xomi + x2li - g54) + rec.d5433 * sin(-xomi + x2li - g54);
                xldot = rec.xni + rec.xfact;
                xnddt = rec.d2201 * cos(x2omi + rec.xli - g22) + rec.d2211 * cos(rec.xli - g22) + ...
                        rec.d3210 * cos(xomi + rec.xli - g32) + rec.d3222 * cos(-xomi + rec.xli - g32) + ...
                        rec.d5220 * cos(xomi + rec.xli - g52) + rec.d5232 * cos(-xomi + rec.xli - g52) + ...
                    2.0 * (rec.d4410 * cos(x2omi + x2li - g44) + ...
                            rec.d4422 * cos(x2li - g44) + rec.d5421 * cos(xomi + x2li - g54) + ...
                            rec.d5433 * cos(-xomi + x2li - g54));
                xnddt = xnddt * xldot;
            end

            %/* ----------------------- integrator ------------------- */
            % sgp4fix move end checks to end of routine
            if (abs(rec.t - rec.atime) >= stepp)
                iretn = 381;
            else % exit here
                ft = rec.t - rec.atime;
                iretn = 0;
            end

            if (iretn == 381)
                rec.xli = rec.xli + xldot * delt + xndt * step2;
                rec.xni = rec.xni + xndt * delt + xnddt * step2;
                rec.atime = rec.atime + delt;
            end
        end

        rec.nm = rec.xni + xndt * ft + xnddt * ft * ft * 0.5;
        xl = rec.xli + xldot * ft + xndt * ft * ft * 0.5;
        if (rec.irez ~= 1)
            rec.mm = xl - 2.0 * rec.nodem + 2.0 * theta;
            rec.dndt = rec.nm - rec.no_unkozai;
        else
            rec.mm = xl - rec.nodem - rec.argpm + theta;
            rec.dndt = rec.nm - rec.no_unkozai;
        end
        rec.nm = rec.no_unkozai + rec.dndt;
    end
end
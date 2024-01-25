function dsinit(tc, xpidot, rec)
    q22 = 1.7891679e-6;
    q31 = 2.1460748e-6;
    q33 = 2.2123015e-7;
    root22 = 1.7891679e-6;
    root44 = 7.3636953e-9;
    root54 = 2.1765803e-9;
    rptim = 4.37526908801129966e-3; % this equates to 7.29211514668855e-5 rad/sec
    root32 = 3.7393792e-7;
    root52 = 1.1428639e-7;
    x2o3 = 2.0 / 3.0;
    znl = 1.5835218e-4;
    zns = 1.19459e-5;

    % sgp4fix identify constants and allow alternate values
    % just xke is used here so pass it in rather than have multiple calls
    % getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

    %/* -------------------- deep space initialization ------------ */
    rec.irez = 0;
    if ((rec.nm < 0.0052359877) && (rec.nm > 0.0034906585))
        rec.irez = 1;
    end

    if ((rec.nm >= 8.26e-3) && (rec.nm <= 9.24e-3) && (rec.em >= 0.5))
        rec.irez = 2;
    end

    %/* ------------------------ do solar terms ------------------- */
    ses = rec.ss1 * zns * rec.ss5;
    sis = rec.ss2 * zns * (rec.sz11 + rec.sz13);
    sls = -zns * rec.ss3 * (rec.sz1 + rec.sz3 - 14.0 - 6.0 * rec.emsq);
    sghs = rec.ss4 * zns * (rec.sz31 + rec.sz33 - 6.0);
    shs = -zns * rec.ss2 * (rec.sz21 + rec.sz23);

    % sgp4fix for 180 deg incl
    if ((rec.inclm < 5.2359877e-2) || (rec.inclm > pi - 5.2359877e-2))
        shs = 0.0;
    end
    if (rec.sinim ~= 0.0)
        shs = shs / rec.sinim;
    end

    sgs = sghs - rec.cosim * shs;

    %/* ------------------------- do lunar terms ------------------ */
    rec.dedt = ses + rec.s1 * znl * rec.s5;
    rec.didt = sis + rec.s2 * znl * (rec.z11 + rec.z13);
    rec.dmdt = sls - znl * rec.s3 * (rec.z1 + rec.z3 - 14.0 - 6.0 * rec.emsq);
    sghl = rec.s4 * znl * (rec.z31 + rec.z33 - 6.0);
    shll = -znl * rec.s2 * (rec.z21 + rec.z23);

    % sgp4fix for 180 deg incl
    if ((rec.inclm < 5.2359877e-2) || (rec.inclm > pi - 5.2359877e-2))
        shll = 0.0;
    end

    rec.domdt = sgs + sghl;
    rec.dnodt = shs;
    
    if (rec.sinim ~= 0.0)
        rec.domdt = rec.domdt - rec.cosim / rec.sinim * shll;
        rec.dnodt = rec.dnodt + shll / rec.sinim;
    end

    %/* ----------- calculate deep space resonance effects -------- */
    rec.dndt = 0.0;
    theta = SGP4.SGP4.fmod(rec.gsto + tc * rptim, SGP4.SGP4.twopi);
    rec.em = rec.em + rec.dedt * rec.t;
    rec.inclm = rec.inclm + rec.didt * rec.t;
    rec.argpm = rec.argpm + rec.domdt * rec.t;
    rec.nodem = rec.nodem + rec.dnodt * rec.t;
    rec.mm = rec.mm + rec.dmdt * rec.t;
    %   sgp4fix for negative inclinations
    %   the following if statement should be commented out
    %if (inclm < 0.0)
    %  {
    %    inclm  = -inclm;
    %    argpm  = argpm - pi;
    %    nodem = nodem + pi;
    %  }

    %/* -------------- initialize the resonance terms ------------- */
    if (rec.irez ~= 0)
        aonv = (rec.nm / rec.xke) ^ x2o3;

        %/* ---------- geopotential resonance for 12 hour orbits ------ */
        if (rec.irez == 2)
            cosisq = rec.cosim * rec.cosim;
            emo = rec.em;
            rec.em = rec.ecco;
            emsqo = rec.emsq;
            rec.emsq = rec.eccsq;
            eoc = rec.em * rec.emsq;
            g201 = -0.306 - (rec.em - 0.64) * 0.440;

            if (rec.em <= 0.65)
                g211 = 3.616 - 13.2470 * rec.em + 16.2900 * rec.emsq;
                g310 = -19.302 + 117.3900 * rec.em - 228.4190 * rec.emsq + 156.5910 * eoc;
                g322 = -18.9068 + 109.7927 * rec.em - 214.6334 * rec.emsq + 146.5816 * eoc;
                g410 = -41.122 + 242.6940 * rec.em - 471.0940 * rec.emsq + 313.9530 * eoc;
                g422 = -146.407 + 841.8800 * rec.em - 1629.014 * rec.emsq + 1083.4350 * eoc;
                g520 = -532.114 + 3017.977 * rec.em - 5740.032 * rec.emsq + 3708.2760 * eoc;
            else
                g211 = -72.099 + 331.819 * rec.em - 508.738 * rec.emsq + 266.724 * eoc;
                g310 = -346.844 + 1582.851 * rec.em - 2415.925 * rec.emsq + 1246.113 * eoc;
                g322 = -342.585 + 1554.908 * rec.em - 2366.899 * rec.emsq + 1215.972 * eoc;
                g410 = -1052.797 + 4758.686 * rec.em - 7193.992 * rec.emsq + 3651.957 * eoc;
                g422 = -3581.690 + 16178.110 * rec.em - 24462.770 * rec.emsq + 12422.520 * eoc;
                if (rec.em > 0.715)
                    g520 = -5149.66 + 29936.92 * rec.em - 54087.36 * rec.emsq + 31324.56 * eoc;
                else
                    g520 = 1464.74 - 4664.75 * rec.em + 3763.64 * rec.emsq;
                end
            end
            
            if (rec.em < 0.7)
                g533 = -919.22770 + 4988.6100 * rec.em - 9064.7700 * rec.emsq + 5542.21  * eoc;
                g521 = -822.71072 + 4568.6173 * rec.em - 8491.4146 * rec.emsq + 5337.524 * eoc;
                g532 = -853.66600 + 4690.2500 * rec.em - 8624.7700 * rec.emsq + 5341.4  * eoc;
            else
                g533 = -37995.780 + 161616.52 * rec.em - 229838.20 * rec.emsq + 109377.94 * eoc;
                g521 = -51752.104 + 218913.95 * rec.em - 309468.16 * rec.emsq + 146349.42 * eoc;
                g532 = -40023.880 + 170470.89 * rec.em - 242699.48 * rec.emsq + 115605.82 * eoc;
            end

            sini2 = rec.sinim * rec.sinim;
            f220 = 0.75 * (1.0 + 2.0 * rec.cosim + cosisq);
            f221 = 1.5 * sini2;
            f321 = 1.875 * rec.sinim  *  (1.0 - 2.0 * rec.cosim - 3.0 * cosisq);
            f322 = -1.875 * rec.sinim  *  (1.0 + 2.0 * rec.cosim - 3.0 * cosisq);
            f441 = 35.0 * sini2 * f220;
            f442 = 39.3750 * sini2 * sini2;
            f522 = 9.84375 * rec.sinim * (sini2 * (1.0 - 2.0 * rec.cosim - 5.0 * cosisq) + ...
                0.33333333 * (-2.0 + 4.0 * rec.cosim + 6.0 * cosisq));
            f523 = rec.sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * rec.cosim + ...
                10.0 * cosisq) + 6.56250012 * (1.0 + 2.0 * rec.cosim - 3.0 * cosisq));
            f542 = 29.53125 * rec.sinim * (2.0 - 8.0 * rec.cosim + cosisq * ...
                (-12.0 + 8.0 * rec.cosim + 10.0 * cosisq));
            f543 = 29.53125 * rec.sinim * (-2.0 - 8.0 * rec.cosim + cosisq * ...
                (12.0 + 8.0 * rec.cosim - 10.0 * cosisq));
            xno2 = rec.nm * rec.nm;
            ainv2 = aonv * aonv;
            temp1 = 3.0 * xno2 * ainv2;
            temp = temp1 * root22;
            rec.d2201 = temp * f220 * g201;
            rec.d2211 = temp * f221 * g211;
            temp1 = temp1 * aonv;
            temp = temp1 * root32;
            rec.d3210 = temp * f321 * g310;
            rec.d3222 = temp * f322 * g322;
            temp1 = temp1 * aonv;
            temp = 2.0 * temp1 * root44;
            rec.d4410 = temp * f441 * g410;
            rec.d4422 = temp * f442 * g422;
            temp1 = temp1 * aonv;
            temp = temp1 * root52;
            rec.d5220 = temp * f522 * g520;
            rec.d5232 = temp * f523 * g532;
            temp = 2.0 * temp1 * root54;
            rec.d5421 = temp * f542 * g521;
            rec.d5433 = temp * f543 * g533;
            rec.xlamo = SGP4.SGP4.fmod(rec.mo + rec.nodeo + rec.nodeo - theta - theta, SGP4.SGP4.twopi);
            rec.xfact = rec.mdot + rec.dmdt + 2.0 * (rec.nodedot + rec.dnodt - rptim) - rec.no_unkozai;
            rec.em = emo;
            rec.emsq = emsqo;
        end

        %/* ---------------- synchronous resonance terms -------------- */
        if (rec.irez == 1)
            g200 = 1.0 + rec.emsq * (-2.5 + 0.8125 * rec.emsq);
            g310 = 1.0 + 2.0 * rec.emsq;
            g300 = 1.0 + rec.emsq * (-6.0 + 6.60937 * rec.emsq);
            f220 = 0.75 * (1.0 + rec.cosim) * (1.0 + rec.cosim);
            f311 = 0.9375 * rec.sinim * rec.sinim * (1.0 + 3.0 * rec.cosim) - 0.75 * (1.0 + rec.cosim);
            f330 = 1.0 + rec.cosim;
            f330 = 1.875 * f330 * f330 * f330;
            rec.del1 = 3.0 * rec.nm * rec.nm * aonv * aonv;
            rec.del2 = 2.0 * rec.del1 * f220 * g200 * q22;
            rec.del3 = 3.0 * rec.del1 * f330 * g300 * q33 * aonv;
            rec.del1 = rec.del1 * f311 * g310 * q31 * aonv;
            rec.xlamo = SGP4.SGP4.fmod(rec.mo + rec.nodeo + rec.argpo - theta, SGP4.SGP4.twopi);
            rec.xfact = rec.mdot + xpidot - rptim + rec.dmdt + rec.domdt + rec.dnodt - rec.no_unkozai;
        end

        %/* ------------ for sgp4, initialize the integrator ---------- */
        rec.xli = rec.xlamo;
        rec.xni = rec.no_unkozai;
        rec.atime = 0.0;
        rec.nm = rec.no_unkozai + rec.dndt;
    end
end
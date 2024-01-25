function [ok] = sgp4_propagator(satrec, tsince, rvhand)
    
    %/* ------------------ set mathematical constants --------------- */
    % sgp4fix divisor for divide by zero check on inclination
    % the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
    % 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
    temp4 = 1.5e-12;
    x2o3 = 2.0 / 3.0;

    % sgp4fix identify constants and allow alternate values
    % getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
    vkmpersec = satrec.radiusearthkm * satrec.xke / 60.0;

    %/* --------------------- clear sgp4 error flag ----------------- */
    satrec.t = tsince;
    satrec.error = 0;

    %/* ------- update for secular gravity and atmospheric drag ----- */
    xmdf = satrec.mo + satrec.mdot * satrec.t;
    argpdf = satrec.argpo + satrec.argpdot * satrec.t;
    nodedf = satrec.nodeo + satrec.nodedot * satrec.t;
    satrec.argpm = argpdf;
    satrec.mm = xmdf;
    t2 = satrec.t * satrec.t;
    satrec.nodem = nodedf + satrec.nodecf * t2;
    tempa = 1.0 - satrec.cc1 * satrec.t;
    tempe = satrec.bstar * satrec.cc4 * satrec.t;
    templ = satrec.t2cof * t2;
    mrt = 0;
    
    if (satrec.isimp ~= 1)
        delomg = satrec.omgcof * satrec.t;
        % sgp4fix use mutliply for speed instead of pow
        delmtemp = 1.0 + satrec.eta * cos(xmdf);
        delm = satrec.xmcof * (delmtemp * delmtemp * delmtemp - satrec.delmo);
        temp = delomg + delm;
        satrec.mm = xmdf + temp;
        satrec.argpm = argpdf - temp;
        t3 = t2 * satrec.t;
        t4 = t3 * satrec.t;
        tempa = tempa - satrec.d2 * t2 - satrec.d3 * t3 - satrec.d4 * t4;
        tempe = tempe + satrec.bstar * satrec.cc5 * (sin(satrec.mm) -satrec.sinmao);
        templ = templ + satrec.t3cof * t3 + t4 * (satrec.t4cof + satrec.t * satrec.t5cof);
    end

    satrec.nm = satrec.no_unkozai;
    satrec.em = satrec.ecco;
    satrec.inclm = satrec.inclo;

    if (satrec.method == 'd')
        tc = satrec.t;
        SGP4.SGP4.dspace(tc,satrec);        
    end

    if (satrec.nm <= 0.0)
        satrec.error = 2;
        ok = false;
    else
        ok = true;
        satrec.am = ((satrec.xke / satrec.nm)^x2o3) * tempa * tempa;
        satrec.nm = satrec.xke / (satrec.am^ 1.5);
        satrec.em = satrec.em - tempe;

        % fix tolerance for error recognition
        % sgp4fix am is fixed from the previous nm check
        if ((satrec.em >= 1.0) || (satrec.em < -0.001))   %/* || (am < 0.95)*/)
            satrec.error = 1;
            % sgp4fix to return if there is an error in eccentricity
            ok = false;
        else
            % sgp4fix fix tolerance to avoid a divide by zero
            if (satrec.em < 1.0e-6)
                satrec.em = 1.0e-6;
            end
            satrec.mm = satrec.mm + satrec.no_unkozai * templ;
            xlm = satrec.mm + satrec.argpm + satrec.nodem;
            satrec.emsq = satrec.em * satrec.em;
        
            satrec.nodem = SGP4.SGP4.fmod(satrec.nodem, SGP4.SGP4.twopi);
            satrec.argpm = SGP4.SGP4.fmod(satrec.argpm, SGP4.SGP4.twopi);
            xlm = SGP4.SGP4.fmod(xlm, SGP4.SGP4.twopi);
            satrec.mm = SGP4.SGP4.fmod(xlm - satrec.argpm - satrec.nodem, SGP4.SGP4.twopi);
        
            % sgp4fix recover singly averaged mean elements
            satrec.am = satrec.am;
            satrec.em = satrec.em;
            satrec.im = satrec.inclm;
            satrec.Om = satrec.nodem;
            satrec.om = satrec.argpm;
            satrec.mm = satrec.mm;
            satrec.nm = satrec.nm;
            
            %/* ----------------- compute extra mean quantities ------------- */
            satrec.sinim = sin(satrec.inclm);
            satrec.cosim = cos(satrec.inclm);
        
            %/* -------------------- add lunar-solar periodics -------------- */
            satrec.ep = satrec.em;
            xincp = satrec.inclm;
            satrec.inclp = satrec.inclm;
            satrec.argpp = satrec.argpm;
            satrec.nodep = satrec.nodem;
            satrec.mp = satrec.mm;
            sinip = satrec.sinim;
            cosip = satrec.cosim;
            
            if (satrec.method == 'd')
                SGP4.SGP4.dpper(satrec.e3, satrec.ee2, satrec.peo, satrec.pgho, ...
                        satrec.pho, satrec.pinco, satrec.plo, satrec.se2, ...
                        satrec.se3, satrec.sgh2, satrec.sgh3, satrec.sgh4, ...
                        satrec.sh2, satrec.sh3, satrec.si2, satrec.si3, ...
                        satrec.sl2, satrec.sl3, satrec.sl4, satrec.t, ...
                        satrec.xgh2, satrec.xgh3, satrec.xgh4, satrec.xh2, ...
                        satrec.xh3, satrec.xi2, satrec.xi3, satrec.xl2, ...
                        satrec.xl3, satrec.xl4, satrec.zmol, satrec.zmos,  ...
                        'n', satrec, satrec.operationmode);
                
                xincp = satrec.inclp;
                if (xincp < 0.0)
                    xincp = -xincp;
                    satrec.nodep = satrec.nodep + pi;
                    satrec.argpp = satrec.argpp - pi;
                end

                if ((satrec.ep < 0.0) || (satrec.ep > 1.0))
                    satrec.error = 3;
                    % sgp4fix add return
                    ok = false;
                end
            end
        
            if (ok)
                %/* -------------------- long period periodics ------------------ */
                if (satrec.method == 'd')
                    sinip = sin(xincp);
                    cosip = cos(xincp);
                    satrec.aycof = -0.5*satrec.j3oj2*sinip;
                    % sgp4fix for divide by zero for xincp = 180 deg
                    if (abs(cosip + 1.0) > 1.5e-12)
                        satrec.xlcof = -0.25 * satrec.j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
                    else
                        satrec.xlcof = -0.25 * satrec.j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
                    end
                end
            
                axnl = satrec.ep * cos(satrec.argpp);
                temp = 1.0 / (satrec.am * (1.0 - satrec.ep * satrec.ep));
                aynl = satrec.ep* sin(satrec.argpp) + temp * satrec.aycof;
                xl = satrec.mp + satrec.argpp + satrec.nodep + temp * satrec.xlcof * axnl;
            
                %/* --------------------- solve kepler's equation --------------- */
                u = SGP4.SGP4.fmod(xl - satrec.nodep, SGP4.SGP4.twopi);
                eo1 = u;
                tem5 = 9999.9;
                ktr = 1;
                sineo1 = 0;
                coseo1 = 0;
            
                %   sgp4fix for kepler iteration
                %   the following iteration needs better limits on corrections
                while ((abs(tem5) >= 1.0e-12) && (ktr <= 10))
                    sineo1 = sin(eo1);
                    coseo1 = cos(eo1);
                    tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
                    tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
                    if (abs(tem5) >= 0.95)
                        if(tem5 > 0.0)
                            tem5 = 0.95;
                        else
                            tem5 = -0.95;
                        end
                    end
                    eo1 = eo1 + tem5;
                    ktr = ktr + 1;
                end
            
                %/* ------------- short period preliminary quantities ----------- */
                ecose = axnl*coseo1 + aynl*sineo1;
                esine = axnl*sineo1 - aynl*coseo1;
                el2 = axnl*axnl + aynl*aynl;
                pl = satrec.am*(1.0 - el2);
            
                if (pl < 0.0)
                    satrec.error = 4;
                    % sgp4fix add return
                    ok = false;
                else
                    rl = satrec.am * (1.0 - ecose);
                    rdotl = sqrt(satrec.am) * esine / rl;
                    rvdotl = sqrt(pl) / rl;
                    betal = sqrt(1.0 - el2);
                    temp = esine / (1.0 + betal);
                    sinu = satrec.am / rl * (sineo1 - aynl - axnl * temp);
                    cosu = satrec.am / rl * (coseo1 - axnl + aynl * temp);
                    su = atan2(sinu, cosu);
                    sin2u = (cosu + cosu) * sinu;
                    cos2u = 1.0 - 2.0 * sinu * sinu;
                    temp = 1.0 / pl;
                    temp1 = 0.5 * satrec.j2 * temp;
                    temp2 = temp1 * temp;
            
                    %/* -------------- update for short period periodics ------------ */
                    if (satrec.method == 'd')
                        cosisq = cosip * cosip;
                        satrec.con41 = 3.0*cosisq - 1.0;
                        satrec.x1mth2 = 1.0 - cosisq;
                        satrec.x7thm1 = 7.0*cosisq - 1.0;
                    end
            
                    mrt = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) + 0.5 * temp1 * satrec.x1mth2 * cos2u;
                    su = su - 0.25 * temp2 * satrec.x7thm1 * sin2u;
                    xnode = satrec.nodep + 1.5 * temp2 * cosip * sin2u;
                    xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
                    mvt = rdotl - satrec.nm * temp1 * satrec.x1mth2 * sin2u / satrec.xke;
                    rvdot = rvdotl + satrec.nm * temp1 * (satrec.x1mth2 * cos2u + 1.5 * satrec.con41) / satrec.xke;
            
                    %/* --------------------- orientation vectors ------------------- */
                    sinsu = sin(su);
                    cossu = cos(su);
                    snod = sin(xnode);
                    cnod = cos(xnode);
                    sini = sin(xinc);
                    cosi = cos(xinc);
                    xmx = -snod * cosi;
                    xmy = cnod * cosi;
                    ux = xmx * sinsu + cnod * cossu;
                    uy = xmy * sinsu + snod * cossu;
                    uz = sini * sinsu;
                    vx = xmx * cossu - cnod * sinsu;
                    vy = xmy * cossu - snod * sinsu;
                    vz = sini * cossu;
            
                    %/* --------- position and velocity (in km and km/sec) ---------- */
                    rvhand.r(1) = (mrt * ux)* satrec.radiusearthkm;
                    rvhand.r(2) = (mrt * uy)* satrec.radiusearthkm;
                    rvhand.r(3) = (mrt * uz)* satrec.radiusearthkm;
                    rvhand.v(1) = (mvt * ux + rvdot * vx) * vkmpersec;
                    rvhand.v(2) = (mvt * uy + rvdot * vy) * vkmpersec;
                    rvhand.v(3) = (mvt * uz + rvdot * vz) * vkmpersec;
                end
            
                % sgp4fix for decaying satellites
                if (mrt < 1.0)
                    satrec.error = 6;
                    ok = false;
                else
                    ok = true;
                end
            end
        end
    end
end  
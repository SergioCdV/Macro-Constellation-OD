function [ok] = sgp4init( opsmode, satrec )
    %/* --------------------- local variables ------------------------ */
    
    %cc1sq, cc2, cc3, coef, coef1, cosio4, eeta, etasq, perige, pinvsq, psisq, qzms24;
    %sfour,tc, temp, temp1, temp2, temp3, tsi, xpidot, xhdot1,qzms2t, ss, x2o3, r, v, delmotemp, qzms2ttemp, qzms24temp;
                   
    epoch = (satrec.jdsatepoch + satrec.jdsatepochF) - 2433281.5;

    %/* ------------------------ initialization --------------------- */
    % sgp4fix divisor for divide by zero check on inclination
    % the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
    % 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
    temp4 = 1.5e-12;

    %/* ----------- set all near earth variables to zero ------------ */
    satrec.isimp = 0;   satrec.method = 'n'; satrec.aycof = 0.0;
    satrec.con41 = 0.0; satrec.cc1 = 0.0; satrec.cc4 = 0.0;
    satrec.cc5 = 0.0; satrec.d2 = 0.0; satrec.d3 = 0.0;
    satrec.d4 = 0.0; satrec.delmo = 0.0; satrec.eta = 0.0;
    satrec.argpdot = 0.0; satrec.omgcof = 0.0; satrec.sinmao = 0.0;
    satrec.t = 0.0; satrec.t2cof = 0.0; satrec.t3cof = 0.0;
    satrec.t4cof = 0.0; satrec.t5cof = 0.0; satrec.x1mth2 = 0.0;
    satrec.x7thm1 = 0.0; satrec.mdot = 0.0; satrec.nodedot = 0.0;
    satrec.xlcof = 0.0; satrec.xmcof = 0.0; satrec.nodecf = 0.0;

    %/* ----------- set all deep space variables to zero ------------ */
    satrec.irez = 0;   satrec.d2201 = 0.0; satrec.d2211 = 0.0;
    satrec.d3210 = 0.0; satrec.d3222 = 0.0; satrec.d4410 = 0.0;
    satrec.d4422 = 0.0; satrec.d5220 = 0.0; satrec.d5232 = 0.0;
    satrec.d5421 = 0.0; satrec.d5433 = 0.0; satrec.dedt = 0.0;
    satrec.del1 = 0.0; satrec.del2 = 0.0; satrec.del3 = 0.0;
    satrec.didt = 0.0; satrec.dmdt = 0.0; satrec.dnodt = 0.0;
    satrec.domdt = 0.0; satrec.e3 = 0.0; satrec.ee2 = 0.0;
    satrec.peo = 0.0; satrec.pgho = 0.0; satrec.pho = 0.0;
    satrec.pinco = 0.0; satrec.plo = 0.0; satrec.se2 = 0.0;
    satrec.se3 = 0.0; satrec.sgh2 = 0.0; satrec.sgh3 = 0.0;
    satrec.sgh4 = 0.0; satrec.sh2 = 0.0; satrec.sh3 = 0.0;
    satrec.si2 = 0.0; satrec.si3 = 0.0; satrec.sl2 = 0.0;
    satrec.sl3 = 0.0; satrec.sl4 = 0.0; satrec.gsto = 0.0;
    satrec.xfact = 0.0; satrec.xgh2 = 0.0; satrec.xgh3 = 0.0;
    satrec.xgh4 = 0.0; satrec.xh2 = 0.0; satrec.xh3 = 0.0;
    satrec.xi2 = 0.0; satrec.xi3 = 0.0; satrec.xl2 = 0.0;
    satrec.xl3 = 0.0; satrec.xl4 = 0.0; satrec.xlamo = 0.0;
    satrec.zmol = 0.0; satrec.zmos = 0.0; satrec.atime = 0.0;
    satrec.xli = 0.0; satrec.xni = 0.0;

    %/* ------------------------ earth constants ----------------------- */
    % sgp4fix identify constants and allow alternate values
    % this is now the only call for the constants
    SGP4.SGP4().getgravconst(satrec.whichconst, satrec);

    %-------------------------------------------------------------------------

    satrec.error = 0;
    satrec.operationmode = opsmode;

    % single averaged mean elements
    satrec.am = 0;
    satrec.em = 0;
    satrec.im = 0;
    satrec.Om = 0;
    satrec.mm = 0;
    satrec.nm = 0.0;

    %/* ------------------------ earth constants ----------------------- */
    % sgp4fix identify constants and allow alternate values no longer needed
    % getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
    ss = 78.0 / satrec.radiusearthkm + 1.0;

    % sgp4fix use multiply for speed instead of pow
    qzms2ttemp = (120.0 - 78.0) / satrec.radiusearthkm;
    qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
    x2o3 = 2.0 / 3.0;

    satrec.init = 'y';
    satrec.t = 0.0;

    % sgp4fix remove satn as it is not needed in initl
    
    SGP4.SGP4.initl(epoch,satrec);
    
    satrec.a = (satrec.no_unkozai * satrec.tumin)^(-2.0 / 3.0);
    satrec.alta = satrec.a * (1.0 + satrec.ecco) - 1.0;
    satrec.altp = satrec.a * (1.0 - satrec.ecco) - 1.0;
    satrec.error = 0;

    % sgp4fix remove this check as it is unnecessary
    % the mrt check in sgp4 handles decaying satellite cases even if the starting
    % condition is below the surface of te earth
    %     if (rp < 1.0)
    %       {
    %         satrec.error = 5;
    %       }

    if ((satrec.omeosq >= 0.0) || (satrec.no_unkozai >= 0.0))
        satrec.isimp = 0;

        if (satrec.rp < (220.0 / satrec.radiusearthkm + 1.0))
            satrec.isimp = 1;
        end

        sfour = ss;
        qzms24 = qzms2t;
        perige = (satrec.rp - 1.0) * satrec.radiusearthkm;

       % /* - for perigees below 156 km, s and qoms2t are altered - */
        if (perige < 156.0)
            sfour = perige - 78.0;
            if (perige < 98.0)
                sfour = 20.0;
            end
            % sgp4fix use multiply for speed instead of pow
            qzms24temp = (120.0 - sfour) / satrec.radiusearthkm;
            qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
            sfour = sfour / satrec.radiusearthkm + 1.0;
        end
        pinvsq = 1.0 / satrec.posq;

        tsi = 1.0 / (satrec.ao - sfour);
        satrec.eta = satrec.ao * satrec.ecco * tsi;
        etasq = satrec.eta * satrec.eta;
        eeta = satrec.ecco * satrec.eta;
        psisq = abs(1.0 - etasq);
        coef = qzms24 * (tsi^ 4.0);
        coef1 = coef / (psisq^3.5);
        cc2 = coef1 * satrec.no_unkozai * (satrec.ao * (1.0 + 1.5 * etasq + eeta * ...
            (4.0 + etasq)) + 0.375 * satrec.j2 * tsi / psisq * satrec.con41 * ...
            (8.0 + 3.0 * etasq * (8.0 + etasq)));
        satrec.cc1 = satrec.bstar * cc2;
        cc3 = 0.0;

        if (satrec.ecco > 1.0e-4)
            cc3 = -2.0 * coef * tsi * satrec.j3oj2 * satrec.no_unkozai * satrec.sinio / satrec.ecco;
        end

        satrec.x1mth2 = 1.0 - satrec.cosio2;
        satrec.cc4 = 2.0* satrec.no_unkozai * coef1 * satrec.ao * satrec.omeosq * ...
            (satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecco * ...
            (0.5 + 2.0 * etasq) - satrec.j2 * tsi / (satrec.ao * psisq) * ...
            (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq * ...
            (1.5 - 0.5 * eeta)) + 0.75 * satrec.x1mth2 * ...
            (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * satrec.argpo)));
        satrec.cc5 = 2.0 * coef1 * satrec.ao * satrec.omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
        cosio4 = satrec.cosio2 * satrec.cosio2;
        temp1 = 1.5 * satrec.j2 * pinvsq * satrec.no_unkozai;
        temp2 = 0.5 * temp1 * satrec.j2 * pinvsq;
        temp3 = -0.46875 * satrec.j4 * pinvsq * pinvsq * satrec.no_unkozai;
        satrec.mdot = satrec.no_unkozai + 0.5 * temp1 * satrec.rteosq * satrec.con41 + 0.0625 * ...
            temp2 * satrec.rteosq * (13.0 - 78.0 * satrec.cosio2 + 137.0 * cosio4);
        satrec.argpdot = -0.5 * temp1 * satrec.con42 + 0.0625 * temp2 * ...
            (7.0 - 114.0 * satrec.cosio2 + 395.0 * cosio4) + ...
            temp3 * (3.0 - 36.0 * satrec.cosio2 + 49.0 * cosio4);
        xhdot1 = -temp1 * satrec.cosio;
        satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * satrec.cosio2) + ...
            2.0 * temp3 * (3.0 - 7.0 * satrec.cosio2)) * satrec.cosio;
        xpidot = satrec.argpdot + satrec.nodedot;
        satrec.omgcof = satrec.bstar * cc3 * cos(satrec.argpo);
        satrec.xmcof = 0.0;

        if (satrec.ecco > 1.0e-4)
            satrec.xmcof = -x2o3 * coef * satrec.bstar / eeta;
        end

        satrec.nodecf = 3.5 * satrec.omeosq * xhdot1 * satrec.cc1;
        satrec.t2cof = 1.5 * satrec.cc1;

        % sgp4fix for divide by zero with xinco = 180 deg
        if (abs(satrec.cosio + 1.0) > 1.5e-12)
            satrec.xlcof = -0.25 * satrec.j3oj2 * satrec.sinio * (3.0 + 5.0 * satrec.cosio) / (1.0 + satrec.cosio);
        else
            satrec.xlcof = -0.25 * satrec.j3oj2 * satrec.sinio * (3.0 + 5.0 * satrec.cosio) / temp4;
        end
        satrec.aycof = -0.5 * satrec.j3oj2 * satrec.sinio;
        % sgp4fix use multiply for speed instead of pow
        delmotemp = 1.0 + satrec.eta * cos(satrec.mo);
        satrec.delmo = delmotemp * delmotemp * delmotemp;
        satrec.sinmao = sin(satrec.mo);
        satrec.x7thm1 = 7.0 * satrec.cosio2 - 1.0;

        %/* --------------- deep space initialization ------------- */
        if ((2 * pi / satrec.no_unkozai) >= 225.0)
            satrec.method = 'd';
            satrec.isimp = 1;
            tc = 0.0;
            satrec.inclm = satrec.inclo;

            SGP4.SGP4.dscom(epoch, satrec.ecco, satrec.argpo, tc, satrec.inclo, satrec.nodeo, satrec.no_unkozai,satrec);                
            
            satrec.ep=satrec.ecco;
            satrec.inclp=satrec.inclo;
            satrec.nodep=satrec.nodeo;
            satrec.argpp=satrec.argpo;
            satrec.mp=satrec.mo;

            
            SGP4.SGP4.dpper(satrec.e3, satrec.ee2, satrec.peo, satrec.pgho, ...
                satrec.pho, satrec.pinco, satrec.plo, satrec.se2, ...
                satrec.se3, satrec.sgh2, satrec.sgh3, satrec.sgh4, ...
                satrec.sh2, satrec.sh3, satrec.si2, satrec.si3, ...
                satrec.sl2, satrec.sl3, satrec.sl4, satrec.t, ...
                satrec.xgh2, satrec.xgh3, satrec.xgh4, satrec.xh2, ...
                satrec.xh3, satrec.xi2, satrec.xi3, satrec.xl2, ...
                satrec.xl3, satrec.xl4, satrec.zmol, satrec.zmos, satrec.init,satrec, ...
                satrec.operationmode);


            satrec.ecco=satrec.ep;
            satrec.inclo=satrec.inclp;
            satrec.nodeo=satrec.nodep;
            satrec.argpo=satrec.argpp;
            satrec.mo=satrec.mp;

            satrec.argpm = 0.0;
            satrec.nodem = 0.0;
            satrec.mm = 0.0;
            
            SGP4.SGP4.dsinit(tc, xpidot, satrec);
        end

        %/* ----------- set variables if not deep space ----------- */
        if (satrec.isimp ~= 1)
            cc1sq = satrec.cc1 * satrec.cc1;
            satrec.d2 = 4.0 * satrec.ao * tsi * cc1sq;
            temp = satrec.d2 * tsi * satrec.cc1 / 3.0;
            satrec.d3 = (17.0 * satrec.ao + sfour) * temp;
            satrec.d4 = 0.5 * temp * satrec.ao * tsi * (221.0 * satrec.ao + 31.0 * sfour) * satrec.cc1;
            satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
            satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 * (12.0 * satrec.d2 + 10.0 * cc1sq));
            satrec.t5cof = 0.2 * (3.0 * satrec.d4 + ...
                12.0 * satrec.cc1 * satrec.d3 + ...
                6.0 * satrec.d2 * satrec.d2 + ...
                15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));
        end
    end
    
    rvhand = SGP4.rvhandle();
    SGP4.SGP4().sgp4_propagator(satrec, 0.0, rvhand);
    delete(rvhand);

    satrec.init = 'n';

    % sgp4fix return boolean. satrec.error contains any error codes
    ok = true;
end  
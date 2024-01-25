
function dpper(e3, ee2, peo, pgho, pho, pinco, plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, ...
               t, xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4, zmol, zmos,init, rec, opsmode)    
   % /* ---------------------- Constants ----------------------------- */
    zns = 1.19459e-5;
    zes = 0.01675;
    znl = 1.5835218e-4;
    zel = 0.05490;

   % /* --------------- Calculate time varying periodics ----------- */
    zm = zmos + zns * t;

    % be sure that the initial call has time set to zero
    if (init == 'y')
        zm = zmos;
    end

    zf = zm + 2.0 * zes * sin(zm);
    sinzf = sin(zf);
    f2 = 0.5 * sinzf * sinzf - 0.25;
    f3 = -0.5 * sinzf * cos(zf);
    ses = se2* f2 + se3 * f3;
    sis = si2 * f2 + si3 * f3;
    sls = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
    sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
    shs = sh2 * f2 + sh3 * f3;
    zm = zmol + znl * t;

    if (init == 'y')
        zm = zmol;
    end

    zf = zm + 2.0 * zel * sin(zm);
    sinzf = sin(zf);
    f2 = 0.5 * sinzf * sinzf - 0.25;
    f3 = -0.5 * sinzf * cos(zf);
    sel = ee2 * f2 + e3 * f3;
    sil = xi2 * f2 + xi3 * f3;
    sll = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
    sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
    shll = xh2 * f2 + xh3 * f3;
    pe = ses + sel;
    pinc = sis + sil;
    pl = sls + sll;
    pgh = sghs + sghl;
    ph = shs + shll;
    
    if (init == 'n')
        pe = pe - peo;
        pinc = pinc - pinco;
        pl = pl - plo;
        pgh = pgh - pgho;
        ph = ph - pho;
        rec.inclp = rec.inclp + pinc;
        rec.ep = rec.ep + pe;
        sinip = sin(rec.inclp);
        cosip = cos(rec.inclp);

       % /* ----------------- apply periodics directly ------------ */
        %  sgp4fix for lyddane choice
        %  strn3 used original inclination - this is technically feasible
        %  gsfc used perturbed inclination - also technically feasible
        %  probably best to readjust the 0.2 limit value and limit discontinuity
        %  0.2 rad = 11.45916 deg
        %  use next line for original strn3 approach and original inclination
        %  if (inclo >= 0.2)
        %  use next line for gsfc version and perturbed inclination
        
        if (rec.inclp >= 0.2)
            ph = ph / sinip;
            pgh = pgh - cosip * ph;
            rec.argpp = rec.argpp + pgh;
            rec.nodep = rec.nodep + ph;
            rec.mp = rec.mp + pl;
        else
           % /* ---- apply periodics with lyddane modification ---- */
            sinop = sin(rec.nodep);
            cosop = cos(rec.nodep);
            alfdp = sinip * sinop;
            betdp = sinip * cosop;
            dalf = ph * cosop + pinc * cosip * sinop;
            dbet = -ph * sinop + pinc * cosip * cosop;
            alfdp = alfdp + dalf;
            betdp = betdp + dbet;
            rec.nodep = SGP4.SGP4.fmod(rec.nodep, SGP4.SGP4.twopi);

            %  sgp4fix for afspc written intrinsic functions
            % nodep used without a trigonometric function ahead
            if ((rec.nodep < 0.0) && (opsmode == 'a'))
                rec.nodep = rec.nodep + SGP4.SGP4.twopi;
            end

            xls = rec.mp + rec.argpp + cosip * rec.nodep;
            dls = pl + pgh - pinc * rec.nodep * sinip;
            xls = xls + dls;
            xls = SGP4.SGP4.fmod(xls,SGP4.SGP4.twopi);
            xnoh = rec.nodep;
            rec.nodep = atan2(alfdp, betdp);

            %  sgp4fix for afspc written intrinsic functions
            % nodep used without a trigonometric function ahead
            if ((rec.nodep < 0.0) && (opsmode == 'a'))
                rec.nodep = rec.nodep + SGP4.SGP4.twopi;
            end

            if (abs(xnoh - rec.nodep) > pi)
                if (rec.nodep < xnoh)
                    rec.nodep = rec.nodep + SGP4.SGP4.twopi;
                else
                    rec.nodep = rec.nodep - SGP4.SGP4.twopi;
                end
            end

            rec.mp = rec.mp + pl;
            rec.argpp = xls - rec.mp - cosip * rec.nodep;
        end
    end
end 
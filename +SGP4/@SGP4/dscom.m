
function dscom (epoch, ep, argpp, tc, inclp, nodep, np, rec)
    %/* -------------------------- constants ------------------------- */
    zes = 0.01675;
    zel = 0.05490;
    c1ss = 2.9864797e-6;
    c1l = 4.7968065e-7;
    zsinis = 0.39785416;
    zcosis = 0.91744867;
    zcosgs = 0.1945905;
    zsings = -0.98088458;
    
    rec.nm = np;
    rec.em = ep;
    rec.snodm = sin(nodep);
    rec.cnodm = cos(nodep);
    rec.sinomm = sin(argpp);
    rec.cosomm = cos(argpp);
    rec.sinim = sin(inclp);
    rec.cosim = cos(inclp);
    rec.emsq = rec.em * rec.em;
    betasq = 1.0 - rec.emsq;
    rec.rtemsq = sqrt(betasq);

    %/* ----------------- initialize lunar solar terms --------------- */
    rec.peo = 0.0;
    rec.pinco = 0.0;
    rec.plo = 0.0;
    rec.pgho = 0.0;
    rec.pho = 0.0;
    rec.day = epoch + 18261.5 + tc / 1440.0;

    xnodce = SGP4.SGP4.fmod(4.5236020 - 9.2422029e-4 * rec.day, SGP4.SGP4.twopi);
    stem = sin(xnodce);
    ctem = cos(xnodce);
    zcosil = 0.91375164 - 0.03568096 * ctem;
    zsinil = sqrt(1.0 - zcosil * zcosil);
    zsinhl = 0.089683511 * stem / zsinil;
    zcoshl = sqrt(1.0 - zsinhl * zsinhl);
    rec.gam = 5.8351514 + 0.0019443680 * rec.day;
    zx = 0.39785416 * stem / zsinil;
    zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
    zx = atan2(zx, zy);
    zx = rec.gam + zx - xnodce;
    zcosgl = cos(zx);
    zsingl = sin(zx);

    %/* ------------------------- do solar terms --------------------- */
    zcosg = zcosgs;
    zsing = zsings;
    zcosi = zcosis;
    zsini = zsinis;
    zcosh = rec.cnodm;
    zsinh = rec.snodm;
    cc = c1ss;
    xnoi = 1.0 / rec.nm;
    
    for lsflg = 1:2
        a1 = zcosg * zcosh + zsing * zcosi * zsinh;
        a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
        a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
        a8 = zsing * zsini;
        a9 = zsing * zsinh + zcosg * zcosi * zcosh;
        a10 = zcosg * zsini;
        a2 = rec.cosim * a7 + rec.sinim * a8;
        a4 = rec.cosim * a9 + rec.sinim * a10;
        a5 = -rec.sinim * a7 + rec.cosim * a8;
        a6 = -rec.sinim * a9 + rec.cosim * a10;

        x1 = a1 * rec.cosomm + a2 * rec.sinomm;
        x2 = a3 * rec.cosomm + a4 * rec.sinomm;
        x3 = -a1 * rec.sinomm + a2 * rec.cosomm;
        x4 = -a3 * rec.sinomm + a4 * rec.cosomm;
        x5 = a5 * rec.sinomm;
        x6 = a6 * rec.sinomm;
        x7 = a5 * rec.cosomm;
        x8 = a6 * rec.cosomm;

        rec.z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
        rec.z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
        rec.z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
        rec.z1 = 3.0 *  (a1 * a1 + a2 * a2) + rec.z31 * rec.emsq;
        rec.z2 = 6.0 *  (a1 * a3 + a2 * a4) + rec.z32 * rec.emsq;
        rec.z3 = 3.0 *  (a3 * a3 + a4 * a4) + rec.z33 * rec.emsq;
        rec.z11 = -6.0 * a1 * a5 + rec.emsq *  (-24.0 * x1 * x7 - 6.0 * x3 * x5);
        rec.z12 = -6.0 *  (a1 * a6 + a3 * a5) + rec.emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
        rec.z13 = -6.0 * a3 * a6 + rec.emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
        rec.z21 = 6.0 * a2 * a5 + rec.emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
        rec.z22 = 6.0 *  (a4 * a5 + a2 * a6) + rec.emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
        rec.z23 = 6.0 * a4 * a6 + rec.emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
        rec.z1 = rec.z1 + rec.z1 + betasq * rec.z31;
        rec.z2 = rec.z2 + rec.z2 + betasq * rec.z32;
        rec.z3 = rec.z3 + rec.z3 + betasq * rec.z33;
        rec.s3 = cc * xnoi;
        rec.s2 = -0.5 * rec.s3 / rec.rtemsq;
        rec.s4 = rec.s3 * rec.rtemsq;
        rec.s1 = -15.0 * rec.em * rec.s4;
        rec.s5 = x1 * x3 + x2 * x4;
        rec.s6 = x2 * x3 + x1 * x4;
        rec.s7 = x2 * x4 - x1 * x3;

        %/* ----------------------- do lunar terms ------------------- */
        if (lsflg == 1)
            rec.ss1 = rec.s1;
            rec.ss2 = rec.s2;
            rec.ss3 = rec.s3;
            rec.ss4 = rec.s4;
            rec.ss5 = rec.s5;
            rec.ss6 = rec.s6;
            rec.ss7 = rec.s7;
            rec.sz1 = rec.z1;
            rec.sz2 = rec.z2;
            rec.sz3 = rec.z3;
            rec.sz11 = rec.z11;
            rec.sz12 = rec.z12;
            rec.sz13 = rec.z13;
            rec.sz21 = rec.z21;
            rec.sz22 = rec.z22;
            rec.sz23 = rec.z23;
            rec.sz31 = rec.z31;
            rec.sz32 = rec.z32;
            rec.sz33 = rec.z33;
            zcosg = zcosgl;
            zsing = zsingl;
            zcosi = zcosil;
            zsini = zsinil;
            zcosh = zcoshl * rec.cnodm + zsinhl * rec.snodm;
            zsinh = rec.snodm * zcoshl - rec.cnodm * zsinhl;
            cc = c1l;
        end
    end

    rec.zmol = SGP4.SGP4.fmod(4.7199672 + 0.22997150  * rec.day - rec.gam, SGP4.SGP4.twopi);
    rec.zmos = SGP4.SGP4.fmod(6.2565837 + 0.017201977 * rec.day, SGP4.SGP4.twopi);

    %/* ------------------------ do solar terms ---------------------- */
    rec.se2 = 2.0 * rec.ss1 * rec.ss6;
    rec.se3 = 2.0 * rec.ss1 * rec.ss7;
    rec.si2 = 2.0 * rec.ss2 * rec.sz12;
    rec.si3 = 2.0 * rec.ss2 * (rec.sz13 - rec.sz11);
    rec.sl2 = -2.0 * rec.ss3 * rec.sz2;
    rec.sl3 = -2.0 * rec.ss3 * (rec.sz3 - rec.sz1);
    rec.sl4 = -2.0 * rec.ss3 * (-21.0 - 9.0 * rec.emsq) * zes;
    rec.sgh2 = 2.0 * rec.ss4 * rec.sz32;
    rec.sgh3 = 2.0 * rec.ss4 * (rec.sz33 - rec.sz31);
    rec.sgh4 = -18.0 * rec.ss4 * zes;
    rec.sh2 = -2.0 * rec.ss2 * rec.sz22;
    rec.sh3 = -2.0 * rec.ss2 * (rec.sz23 - rec.sz21);

    %/* ------------------------ do lunar terms ---------------------- */
    rec.ee2 = 2.0 * rec.s1 * rec.s6;
    rec.e3 = 2.0 * rec.s1 * rec.s7;
    rec.xi2 = 2.0 * rec.s2 * rec.z12;
    rec.xi3 = 2.0 * rec.s2 * (rec.z13 - rec.z11);
    rec.xl2 = -2.0 * rec.s3 * rec.z2;
    rec.xl3 = -2.0 * rec.s3 * (rec.z3 - rec.z1);
    rec.xl4 = -2.0 * rec.s3 * (-21.0 - 9.0 * rec.emsq) * zel;
    rec.xgh2 = 2.0 * rec.s4 * rec.z32;
    rec.xgh3 = 2.0 * rec.s4 * (rec.z33 - rec.z31);
    rec.xgh4 = -18.0 * rec.s4 * zel;
    rec.xh2 = -2.0 * rec.s2 * rec.z22;
    rec.xh3 = -2.0 * rec.s2 * (rec.z23 - rec.z21);
end 
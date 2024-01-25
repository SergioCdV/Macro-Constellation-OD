function initl(epoch, rec)    
    %/* ----------------------- earth constants ---------------------- */
    % sgp4fix identify constants and allow alternate values
    % only xke and j2 are used here so pass them in directly
    % getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
    x2o3 = 2.0 / 3.0;

    %/* ------------- calculate auxillary epoch quantities ---------- */
    rec.eccsq = rec.ecco * rec.ecco;
    rec.omeosq = 1.0 - rec.eccsq;
    rec.rteosq = sqrt(rec.omeosq);
    rec.cosio = cos(rec.inclo);
    rec.cosio2 = rec.cosio * rec.cosio;

    %/* ------------------ un-kozai the mean motion ----------------- */
    ak = (rec.xke / rec.no_kozai) ^ x2o3;
    d1 = 0.75 * rec.j2 * (3.0 * rec.cosio2 - 1.0) / (rec.rteosq * rec.omeosq);
    del = d1 / (ak * ak);
    adel = ak * (1.0 - del * del - del * ...
        (1.0 / 3.0 + 134.0 * del * del / 81.0));
    del = d1 / (adel * adel);
    rec.no_unkozai = rec.no_kozai / (1.0 + del);

    rec.ao = (rec.xke / (rec.no_unkozai))^x2o3;
    rec.sinio = sin(rec.inclo);
    po = rec.ao * rec.omeosq;
    rec.con42 = 1.0 - 5.0 * rec.cosio2;
    rec.con41 = -rec.con42 - rec.cosio2 - rec.cosio2;
    rec.ainv = 1.0 / rec.ao;
    rec.posq = po * po;
    rec.rp = rec.ao * (1.0 - rec.ecco);
    rec.method = 'n';

    % find greenwich location at epoch
    rec.gsto = SGP4.SGP4.gstime(epoch + 2433281.5);
end 
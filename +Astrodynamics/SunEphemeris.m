
function rs = SunEphemeris(JD)
    % Orbital elements
    T = (JD-2451545.0)/36525.0;
    M = deg2rad(357.5277233 + 35999.05034 * T);
    Omega = deg2rad(280.460 + 36000.771 * T);

    % Geocentric coordinates 
    eps = deg2rad(23.43929111 - 0.0130042 * T);
    lambda = Omega + deg2rad(6892/3600)*sin(M)+deg2rad(72/3600)*sin(2*M);
    r = (149.619-2.499*cos(M)-0.021*cos(2*M)) * 1e9;
    rs = r .* [cos(lambda); sin(lambda) .* cos(eps); sin(lambda) .* sin(eps)];
end
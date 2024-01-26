


% J2 osculating dynamics 
function [dS] = APSO_dynamics(mu, J2, Re, s, InitialEpoch, tspan)
   r = s(1:3);      % ECI position coordinates
   v = s(4:6);      % ECI velocity coordinates 

   % Dynamics
   R = norm(r);
   gamma = -mu*r/R^3.*[1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-1); ...
                       1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-1); ...
                       1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-3)];
   dS = [v; gamma];
end
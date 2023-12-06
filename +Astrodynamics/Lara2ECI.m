

function [s] = Lara2ECI(H, x, direction)
    if (direction)
        % Polar variables 
        psi = x(1,1); 
        chi = x(2,1); 
        xi = x(3,1);
        r = x(4,1);
        R = x(5,1);
        Theta = x(6,1);

        c = H/Theta;
        
        t = 1-chi^2/(1+c);
        tau = 1-xi^2/(1+c);
        q = chi * xi/(1+c);

        s(1,1) = r * (t*cos(psi)+q*sin(psi));
        s(2,1) = r * (t*sin(psi)-q*cos(psi));
        s(3,1) = r * chi;
        s(4,1) = s(1,1) * R/r - Theta/r * (q*cos(psi)+tau*sin(psi));
        s(5,1) = s(2,1) * R/r - Theta/r * (q*sin(psi)-tau*cos(psi));
        s(6,1) = s(3,1) * R/r + Theta/r * xi;

    else
        s(4,1) = norm(x(1:3,1));
        s(5,1) = dot(x(1:3,1),x(4:6,1))/s(1,1);
        s(6,1) = norm(cross(x(1:3,1),x(4:6,1)));
        s(3,1) = (s(4,1)*x(6,1)-x(3,1)*s(5,1))/s(6,1);
        s(2,1) = x(3,1)/s(1,1);
        
        c = sqrt(1-s(2,1)^2-s(3,1)^2);
        t  = 1-s(2,1)/(1+c);
        q = s(3,1)*s(2,1)/(1+c);

        s(1,1) = atan2((x(1,1)*q+x(2,1)*t)/(s(1,1)*(t^2+q^2)), (x(1,1)*t-x(2,1)*q)/(s(1,1)*(t^2+q^2)));
    end
end
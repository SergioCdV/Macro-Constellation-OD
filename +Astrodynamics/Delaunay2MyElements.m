
function [X] = Delaunay2MyElements(x, direction)
    if (direction)
        % Delaunay elements of the particle
        M = x(1,1);
        omega = x(2,1);
        Omega = x(3,1);
        L = x(4,1);
        G = x(5,1);
        H = x(6,1);

        cos_i = H/G;

        q(1,1) = sin(acos(cos_i)/2) * cos((Omega-omega)/2);
        q(2,1) = sin(acos(cos_i)/2) * sin((Omega-omega)/2);
        q(3,1) = cos(acos(cos_i)/2) * sin((Omega+omega)/2);
        q(4,1) = cos(acos(cos_i)/2) * cos((Omega+omega)/2);

        % Assemble the set
        X = [q; L; G; H; M];  

    else
        % Elements of the particle 
        qp = x(1:4,1);
        L = x(5,1);
        G = x(6,1);
        H = x(7,1);
    
        % Compute the RAAN and AoP from qp 
        diff = atan2(qp(2,1), qp(1,1));
        plus = atan2(qp(3,1), qp(4,1));
        Omega = (plus+diff);
        omega = 2 * plus - Omega;
    
        % Assemble the set
        X = [0 omega Omega L G H].';
    end
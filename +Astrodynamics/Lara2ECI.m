

function [s] = Lara2ECI(H, x, direction)
    if (direction)
        % Polar variables 
        psi = x(1,:); 
        chi = x(2,:); 
        xi = x(3,:);
        r = x(4,:);
        R = x(5,:);
        Theta = x(6,:);

        c = H ./ Theta;
        
        b = 1 - chi.^2 ./ (1+c);
        tau = 1 - xi.^2. / (1+c);
        q = chi .* xi ./ (1+c);

        s(1,:) = r .* (b.*cos(psi)+q.*sin(psi));
        s(2,:) = r .* (b.*sin(psi)-q.*cos(psi));
        s(3,:) = r .* chi;
        s(4,:) = s(1,:) .* R./r - Theta./r .* (q.*cos(psi)+tau.*sin(psi));
        s(5,:) = s(2,:) .* R./r - Theta./r .* (q.*sin(psi)-tau.*cos(psi));
        s(6,:) = s(3,:) .* R./r + Theta./r .* xi;

    else
        s(4,:) = sqrt( dot(x(1:3,:), x(1:3,:), 1) );
        s(5,:) = dot(x(1:3,:), x(4:6,:), 1) ./ s(4,:);
        h = cross( x(1:3,:), x(4:6,:) );
        s(6,:) = sqrt( dot(h,h,1) );
        s(3,:) = (s(4,:) .* x(6,:) - x(3,:) .* s(5,:)) ./ s(6,:);
        s(2,:) = x(3,:) ./ s(4,:);
        
        c = h(3,:) ./ s(6,:);
        t = 1 - s(2,:).^2 ./ (1+c); 
        q = s(3,:) .* s(2,:) ./ (1+c);
        
        den = s(4,:) .* (t.^2 + q.^2);
        s(1,:) = atan2( (x(1,:) .* q + x(2,:) .* t) ./ den, (x(1,:) .* t - x(2,:) .* q) ./ den);
    end
end
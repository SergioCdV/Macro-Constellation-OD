
function [Q] = left_isoclinic(q)
    Q(1,:) = [q(1) -q(2:4).'];
    Q(2,:) = [q(2) q(1) -q(4) q(3)];
    Q(3,:) = [q(3) q(4) q(1) -q(2)];
    Q(4,:) = [q(4) -q(3) q(2) q(1)];
end
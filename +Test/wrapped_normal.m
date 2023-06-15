%% Combinatorial constellation macro-determination 
% Date: 27/10/2022
% Author: Sergio Cuevas del Valle

%% Input data 
mu = deg2rad(0);                % Mean anomaly 
sigma = [8e-3, 0.02, 0.05 0.1];   % Uncertainty parameter 
M = linspace(0,2*pi,1e3);       % Anomaly space 
eps = 1e-7;                     % Tolerance error

MTT = Filters.MTT_filter(zeros(6,1), 1, 5e1, 1, 1, 0);

%% Compute the wrapped normal 
f = zeros(size(sigma,2), size(M,2));
for i = 1:size(sigma,2)
    f(i,:) = MTT.wrapped_normal(eps, M.', mu, sigma(i)) + MTT.wrapped_normal(eps, M.', pi/2, sigma(i));
end

% Results
figure 
hold on
view(3)
plot3(cos(M), sin(M), zeros(1,length(M)))
for i = 1:size(sigma,2)
    plot3(cos(M), sin(M), f(i,:))
end
hold off
zlabel('$f(x)$')
grid on;
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));

%% %% Compute the von Mises 
mu = [1 1 1] / sqrt(3);     % Mean direction
kappa = [1e-3 1e-1 10 50];    % Distribution parameter 
N = 300;                    % Number of samples 
p = size(mu,2);

figure 
view(3)
hold on
quiver3(zeros(1,1), zeros(1,1), zeros(1,1), mu(:,1), mu(:,2), mu(:,3), 'k')
for j = 1:size(kappa,2)
    tmpMu = [1 zeros(1,p-1)];
    t = randVMFMeanDir(N, kappa(j), p);
    
    % RandSphere = QuaternionAlgebra.UniformSphere(N).';
    randNorm = normrnd(zeros(N,p-1), 1, [N, p-1]);
    RandSphere = zeros(N,p-1);
    for r=1:N
        RandSphere(r,:) = randNorm(r,:)./norm(randNorm(r,:));
    end
    RandVMF = repmat(t, [1 p]).*repmat(tmpMu, [N 1]) + repmat((1-t.^2).^(1/2), [1 p]).*[zeros(N,1) RandSphere ];
    
    % Rotate the distribution to the right direction
    Otho = null(mu);
    Rot = [mu' Otho];
    RandVMF = (Rot*RandVMF')';
    scatter3(RandVMF(:,1), RandVMF(:,2), RandVMF(:,3), 'filled');
end
legend('$\mu$', '$\kappa = 0.001$', '$\kappa = 0.1$', '$\kappa = 10$', '$\kappa = 50$', 'AutoUpdate', 'off')
m = 100;
[aa, bb, cc] = sphere(m);
h = surf(aa, bb, cc);
set(h, 'FaceColor',[0 0 1], 'FaceAlpha',0.1,'FaceLighting','gouraud','EdgeColor','none')
axis('equal')
hold off
grid on;
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));

%% Auxiliary functions 



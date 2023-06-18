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
scatter3(mu(:,1), mu(:,2), mu(:,3), 'filled', 'k')
legend('$\kappa = 0.001$', '$\kappa = 0.1$', '$\kappa = 10$', '$\kappa = 50$', '$\mu$', 'AutoUpdate', 'off')
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
function [t] = randVMFMeanDir(N, k, p)
% This function generate random samples from the tangent direction of von
% Mises-Fisher distribution using rejection sampling. The density of is
% described in VMFMeanDirDensity function. See the
% SphereDistributionsRand.pdf file for detail description and formulas.
%
% Usage:
%   [t] = randVMFMeanDir(N, k, p);
%
% Inputs:
%   N: The number of samples one wants to generate.
%
%   k: The kappa parameter of the VMF distribution.
%
%   p: The dimension of the VMF distribution.
%
% Outputs:
%   t : A N x 1 vector which are the random samples from the VMF's tangent
%   distribution.
%
% Function is written by Yu-Hui Chen, University of Michigan
% Contact E-mail: yuhuic@umich.edu
%
min_thresh = 1/(5*N);
xx = -1:0.000001:1;
yy = VMFMeanDirDensity(xx, k, p);
cumyy = cumsum(yy)*(xx(2)-xx(1));
leftBound = xx(find(cumyy>min_thresh,1));
%%% Fin the left bound
xx = linspace(leftBound, 1, 1000);
yy = VMFMeanDirDensity(xx, k, p);
M = max(yy);
t = zeros(N,1);
for i=1:N
    while(1)
        x = rand*(1-leftBound)+leftBound;
        h = VMFMeanDirDensity(x, k, p);
        draw = rand*M;
        if(draw<=h)
            break;
        end
    end
    t(i) = x;
end
end

function [y]=VMFMeanDirDensity(x, k, p)
% This is the tangent direction density of VMF distribution. See the
% SphereDistributionsRand.pdf file for detail description and formulas.
%
% Usage:
%   [y]=VMFMeanDirDensity(x, k, p);
%
% Inputs:
%   x: The tangent direction value. should be in [-1 1].
%
%   k: The kappa parameter of the VMF distribution.
%
%   p: The dimension of the VMF distribution.
%
% Outputs:
%   y : The density value of the VMF tangent density.
%
% Function is written by Yu-Hui Chen, University of Michigan
% Contact E-mail: yuhuic@umich.edu
%
if(any((x<-1) | (x>1)))
    error('Input of x should be within -1~1');
end
Coeff = (k/2)^(p/2-1) * (gamma((p-1)/2)*gamma(1/2)*besseli(p/2-1,k))^(-1);
y = Coeff * exp(k*x).*(1-x.^2).^((p-3)/2);
end

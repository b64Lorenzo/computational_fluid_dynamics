%% =========================================================================
% EXERCISE 2(b)
%
% METHOD A + IMPLICIT EULER
%
% Stability Analysis
%
% h+/h- = 0.95
%
% Tasks:
%
%   1. Build the semi-discrete operator
%   2. Compute the eigenvalue spectrum
%   3. Verify implicit Euler amplification factors
%   4. Verify numerically with transient simulations
%
%% =========================================================================

clear
clc
close all

%% =========================================================================
% Logging
%% =========================================================================

problem_name = 'exercise2_methodA_implicit';

fid = logger_init(problem_name);

logger(fid,'INFO','Starting Exercise 2(b)');

%% =========================================================================
% Directories
%% =========================================================================

results_dir = fullfile('results',problem_name);

figures_dir = fullfile(results_dir,'figures');
data_dir    = fullfile(results_dir,'data');

if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

if ~exist(data_dir,'dir')
    mkdir(data_dir);
end

%% =========================================================================
% Parameters
%% =========================================================================

L = 1.0;
U = 1.0;
k = 1e-2;

N = 101;

stretch_ratio = 0.95;

t_end = 0.01;

%% =========================================================================
% Mesh
%% =========================================================================

r = stretch_ratio;

dx0 = L*(1-r)/(1-r^N);

dx = dx0*r.^(0:N-1);

x = [0; cumsum(dx(:))];

mesh.x = x;

%% =========================================================================
% Mesh figure
%% =========================================================================

fmesh = figure('Color','w');

plot(x,...
     zeros(size(x)),...
     'o')

xlabel('x')
title(sprintf('Method A Mesh (h+/h- = %.2f)',stretch_ratio))

grid on
box on

saveas(fmesh,...
    fullfile(figures_dir,...
    'mesh.png'));

%% =========================================================================
% Boundary Conditions
%% =========================================================================

bc.W.type  = 'dirichlet';
bc.W.value = 0.0;

bc.E.type  = 'dirichlet';
bc.E.value = 1.0;

%% =========================================================================
% Initial Condition
%% =========================================================================

bc.IC = sin(pi*x);

%% =========================================================================
% Semi-discrete operator
%% =========================================================================

n_internal = length(x)-2;

A = zeros(n_internal);

for I = 1:n_internal

    i = I + 1;

    [aW,aP,aE,~] = ...
        convection_methodA_1d( ...
        mesh,U,k,i);

    %
    % Operator consistent with transient formulation
    %

    A(I,I) = -aP;

    if I > 1
        A(I,I-1) = -aW;
    end

    if I < n_internal
        A(I,I+1) = -aE;
    end

end

%% =========================================================================
% Eigenvalues
%% =========================================================================

lambda = eig(A);

max_real = max(real(lambda));
min_real = min(real(lambda));

fprintf('\n');
fprintf('====================================================\n');
fprintf('EIGENVALUE ANALYSIS\n');
fprintf('====================================================\n');

fprintf('min Re(lambda) = %.6e\n',min_real);
fprintf('max Re(lambda) = %.6e\n',max_real);

fprintf('====================================================\n');

logger(fid,'INFO',...
       'min Re(lambda) = %.6e',min_real);

logger(fid,'INFO',...
       'max Re(lambda) = %.6e',max_real);

%% =========================================================================
% Figure 1
% Eigenvalue Spectrum
%% =========================================================================

f1 = figure('Color','w');

plot(real(lambda),...
     imag(lambda),...
     'bo')

hold on

xline(0,...
      'k--',...
      'Re(\lambda)=0')

xlabel('Re(\lambda)')
ylabel('Im(\lambda)')

title('Method A Eigenvalue Spectrum')

grid on
box on

saveas(f1,...
    fullfile(figures_dir,...
    'eigenvalue_spectrum.png'));

%% =========================================================================
% Implicit Euler Amplification Factor
%
% g = 1/(1-dt*lambda)
%
%% =========================================================================

dt_list = logspace(-8,1,200);

rho = zeros(size(dt_list));

for n = 1:length(dt_list)

    dt = dt_list(n);

    g = 1./(1 - dt*lambda);

    rho(n) = max(abs(g));

end

%% =========================================================================
% Figure 2
% Spectral Radius
%% =========================================================================

f2 = figure('Color','w');

loglog(dt_list,...
       rho,...
       'LineWidth',2)

hold on

yline(1,...
      'k--',...
      '|\it g|=1')

xlabel('\Delta t')
ylabel('\rho')

title('Implicit Euler Amplification Factor')

grid on
box on

saveas(f2,...
    fullfile(figures_dir,...
    'amplification_factor.png'));

%% =========================================================================
% Numerical Verification
%% =========================================================================

dt_test = [ ...
    1e-8 ...
    1e-6 ...
    1e-4 ...
    1e-3 ...
    1e-2 ...
    1e-1 ...
    1 ];

energy_ratio = zeros(size(dt_test));

for iCase = 1:length(dt_test)

    dt = dt_test(iCase);

    try

        [phi,time] = ...
            unsteady_convection_diffusion_1d( ...
            mesh,...
            U,...
            k,...
            bc,...
            @convection_methodA_1d,...
            dt,...
            t_end,...
            'implicit_euler',...
            fid);

        energy_ratio(iCase) = ...
            norm(phi(:,end),2) ...
           / norm(phi(:,1),2);

    catch

        energy_ratio(iCase) = Inf;

    end

end

%% =========================================================================
% Figure 3
% Numerical Verification
%% =========================================================================

f3 = figure('Color','w');

semilogx(dt_test,...
         energy_ratio,...
         'o-',...
         'LineWidth',2,...
         'MarkerSize',8)

hold on

yline(1,...
      'k--')

xlabel('\Delta t')

ylabel('||\phi(T)||_2 / ||\phi(0)||_2')

title('Implicit Euler Verification')

grid on
box on

saveas(f3,...
    fullfile(figures_dir,...
    'implicit_verification.png'));

%% =========================================================================
% Save Results
%% =========================================================================

save(fullfile(data_dir,...
     'exercise2_methodA_implicit_results.mat'),...
     'lambda',...
     'dt_list',...
     'rho',...
     'dt_test',...
     'energy_ratio');

%% =========================================================================
% Summary
%% =========================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('FINAL RESULTS\n');
fprintf('====================================================\n');

fprintf('min Re(lambda) = %.6e\n',min_real);
fprintf('max Re(lambda) = %.6e\n',max_real);

fprintf('\n');

fprintf('max spectral radius = %.6e\n',max(rho));

fprintf('\n');

if max_real < 0

    fprintf('All eigenvalues lie in the left half-plane.\n');
    fprintf('Implicit Euler is unconditionally stable.\n');

else

    fprintf('Positive real eigenvalues detected.\n');
    fprintf('Unconditional stability not guaranteed.\n');

end

fprintf('====================================================\n');

logger(fid,'INFO','Exercise completed');

logger_close(fid,problem_name);
%% =========================================================================
% EXERCISE 2(a)
%
% METHOD A + EXPLICIT EULER
%
% Stability Analysis
%
% h+/h- = 0.95
%
% Performs:
%
% 1. Positive coefficient analysis
% 2. Eigenvalue analysis
% 3. Explicit Euler stability estimate
% 4. Numerical verification with transient solver
% 5. Comparison against uniform mesh
%
%% =========================================================================

clear
clc
close all

%% =========================================================================
% Logging
%% =========================================================================

problem_name = 'exercise2_methodA_explicit';

fid = logger_init(problem_name);

logger(fid,'INFO','Starting Exercise 2(a)');

%% =========================================================================
% Output folders
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

t_end = 1e-3;

%% =========================================================================
% Boundary Conditions
%% =========================================================================

bc.W.type  = 'dirichlet';
bc.W.value = 0.0;

bc.E.type  = 'dirichlet';
bc.E.value = 1.0;

%% =========================================================================
% Uniform Mesh
%% =========================================================================

x_uni = linspace(0,L,N)';

mesh_uni.x = x_uni;

%% =========================================================================
% Stretched Mesh
%% =========================================================================

r = stretch_ratio;

h0 = L*(1-r)/(1-r^(N-1));

dx = h0*r.^(0:N-2);

x_str = [0 cumsum(dx)]';

mesh_str.x = x_str;

%% =========================================================================
% Mesh figure
%% =========================================================================

fmesh = figure('Color','w');

subplot(2,1,1)

plot(x_uni,...
     zeros(size(x_uni)),...
     'o')

title('Uniform Mesh')

grid on
box on

subplot(2,1,2)

plot(x_str,...
     zeros(size(x_str)),...
     'o')

title(sprintf('Stretched Mesh (h+/h- = %.2f)',...
      stretch_ratio))

grid on
box on

saveas(fmesh,...
    fullfile(figures_dir,...
    'mesh_comparison.png'));

%% =========================================================================
% Positive coefficient analysis
%% =========================================================================

n_internal = N-2;

aW_all = zeros(n_internal,1);
aP_all = zeros(n_internal,1);
aE_all = zeros(n_internal,1);

dt_pos = inf;

for I = 1:n_internal

    i = I+1;

    [aW,aP,aE,~] = ...
        convection_methodA_1d( ...
        mesh_str,U,k,i);

    aW_all(I)=aW;
    aP_all(I)=aP;
    aE_all(I)=aE;

    if aP < 0
        dt_pos = min(dt_pos,-1/aP);
    end

end

h_min = min(diff(x_str));

eta_pos = U*dt_pos/h_min;

%% =========================================================================
% Report positivity results
%% =========================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('POSITIVE COEFFICIENT ANALYSIS\n');
fprintf('====================================================\n');

fprintf('all(aW>=0) = %d\n',all(aW_all>=0));
fprintf('all(aE>=0) = %d\n',all(aE_all>=0));

fprintf('dt_pos     = %.6e\n',dt_pos);
fprintf('eta_pos    = %.6e\n',eta_pos);

fprintf('====================================================\n');

%% =========================================================================
% Build matrices
%% =========================================================================

A_uni = zeros(n_internal);

for I = 1:n_internal

    i = I+1;

    [aW,aP,aE,~] = ...
        convection_methodA_1d( ...
        mesh_uni,U,k,i);

    %
    % Consistent with transient formulation
    %

    A_uni(I,I) = -aP;

    if I > 1
        A_uni(I,I-1) = -aW;
    end

    if I < n_internal
        A_uni(I,I+1) = -aE;
    end

end

A_str = zeros(n_internal);

for I = 1:n_internal

    i = I+1;

    [aW,aP,aE,~] = ...
        convection_methodA_1d( ...
        mesh_str,U,k,i);

    A_str(I,I) = -aP;

    if I > 1
        A_str(I,I-1) = -aW;
    end

    if I < n_internal
        A_str(I,I+1) = -aE;
    end

end

%% =========================================================================
% Eigenvalues
%% =========================================================================

lambda_uni = eig(A_uni);
lambda_str = eig(A_str);

%% =========================================================================
% Stability limits
%% =========================================================================

dt_candidates = [];

for j=1:length(lambda_uni)

    lam = lambda_uni(j);

    if real(lam)<0

        dt_tmp = ...
            -2*real(lam)/(abs(lam)^2);

        if dt_tmp>0
            dt_candidates(end+1)=dt_tmp;
        end

    end

end

dt_max_uni = min(dt_candidates);

dt_candidates = [];

for j=1:length(lambda_str)

    lam = lambda_str(j);

    if real(lam)<0

        dt_tmp = ...
            -2*real(lam)/(abs(lam)^2);

        if dt_tmp>0
            dt_candidates(end+1)=dt_tmp;
        end

    end

end

dt_max_str = min(dt_candidates);

%% =========================================================================
% Analytical results
%% =========================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('EIGENVALUE STABILITY ANALYSIS\n');
fprintf('====================================================\n');

fprintf('Uniform grid\n');
fprintf('dt_max = %.6e\n',dt_max_uni);

fprintf('\n');

fprintf('Stretched grid\n');
fprintf('dt_max = %.6e\n',dt_max_str);

fprintf('\n');

fprintf('dt_stretch/dt_uniform = %.6e\n',...
        dt_max_str/dt_max_uni);

fprintf('====================================================\n');

%% =========================================================================
% Figure 1
% Eigenvalue spectrum
%% =========================================================================

f1 = figure('Color','w');

plot(real(lambda_uni),...
     imag(lambda_uni),...
     'bo')

hold on

plot(real(lambda_str),...
     imag(lambda_str),...
     'rx')

xline(0,'k--')

xlabel('Re(\lambda)')
ylabel('Im(\lambda)')

title('Method A Eigenvalue Spectrum')

legend('Uniform',...
       'Stretched',...
       'Location','best')

grid on
box on

saveas(f1,...
    fullfile(figures_dir,...
    'eigenvalue_spectrum.png'));

%% =========================================================================
% Figure 2
% Explicit Euler Stability Disk
%% =========================================================================

theta = linspace(0,2*pi,500);

f2 = figure('Color','w');

plot(-1+cos(theta),...
      sin(theta),...
      'k--',...
      'LineWidth',2)

hold on

plot(real(dt_max_uni*lambda_uni),...
     imag(dt_max_uni*lambda_uni),...
     'bo')

plot(real(dt_max_str*lambda_str),...
     imag(dt_max_str*lambda_str),...
     'rx')

axis equal

xlabel('Re(\Delta t \lambda)')
ylabel('Im(\Delta t \lambda)')

title('Explicit Euler Stability Region')

legend('Boundary',...
       'Uniform',...
       'Stretched',...
       'Location','best')

grid on
box on

saveas(f2,...
    fullfile(figures_dir,...
    'explicit_euler_region.png'));

%% =========================================================================
% Numerical Verification
%% =========================================================================

bc.IC = sin(pi*x_str);

dt_test = [ ...
    0.5*dt_max_str ...
    dt_max_str ...
    2.0*dt_max_str ];

final_energy = zeros(size(dt_test));

for m = 1:length(dt_test)

    dt = dt_test(m);

    try

        [phi,time] = ...
            unsteady_convection_diffusion_1d( ...
            mesh_str,...
            U,...
            k,...
            bc,...
            @convection_methodA_1d,...
            dt,...
            t_end,...
            'explicit_euler',...
            fid);

        final_energy(m) = ...
            norm(phi(:,end),2) ...
          / norm(phi(:,1),2);

    catch

        final_energy(m)=Inf;

    end

end

%% =========================================================================
% Figure 3
% Verification
%% =========================================================================

f3 = figure('Color','w');

semilogx(dt_test,...
         final_energy,...
         'o-',...
         'LineWidth',2)

hold on

xline(dt_max_str,...
      'r--',...
      'Predicted dt_{max}')

yline(1,...
      'k--')

xlabel('\Delta t')
ylabel('Energy Amplification')

title('Explicit Euler Verification')

grid on
box on

saveas(f3,...
    fullfile(figures_dir,...
    'explicit_verification.png'));

%% =========================================================================
% Save
%% =========================================================================

save(fullfile(data_dir,...
    'exercise2_results.mat'));

%% =========================================================================
% Final Summary
%% =========================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('FINAL CONCLUSIONS\n');
fprintf('====================================================\n');

fprintf('Positivity dt      = %.6e\n',dt_pos);
fprintf('Eigenvalue dt_max  = %.6e\n',dt_max_str);

fprintf('\n');

fprintf('Uniform dt_max     = %.6e\n',dt_max_uni);
fprintf('Stretched dt_max   = %.6e\n',dt_max_str);

fprintf('\n');

fprintf('Ratio              = %.6e\n',...
        dt_max_str/dt_max_uni);

fprintf('====================================================\n');

logger(fid,'INFO','Exercise completed');

logger_close(fid,problem_name);
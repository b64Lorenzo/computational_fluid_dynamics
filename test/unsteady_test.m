%% =========================================================================
% CFD - 1D TRANSIENT -> STEADY VALIDATION
%
% Verify:
%
%   1) phi(t) -> phi_steady
%
%   2) Residual -> 0
%
%% =========================================================================

clear
clc
close all

%% Logging

problem_name = 'transient_to_steady_validation';

fid = logger_init(problem_name);

logger(fid,'INFO','Starting simulation');

%% =========================================================================
% Output directories
%% =========================================================================

results_dir = fullfile('results',problem_name);

figures_dir = fullfile(results_dir,'figures');

data_dir = fullfile(results_dir,'data');

if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

if ~exist(data_dir,'dir')
    mkdir(data_dir);
end

%% =========================================================================
% Parameters
%% =========================================================================

L = 10.0;
U = 50.0;
k = 1e-3;

N = 100;

dt = 1e-4;
t_end = 0.5;

%% Boundary Conditions

bc.W.type  = 'dirichlet';
bc.W.value = 0.0;

bc.E.type  = 'dirichlet';
bc.E.value = 20.0;

%% Mesh

mesh = create_uniform_mesh_1d(L,N);

x = mesh.x(:);

%% Initial Condition

bc.IC = zeros(size(x));

%% =========================================================================
% Steady Solution
%% =========================================================================

logger(fid,'INFO','Computing steady solution');

[phi_ss,A_ss,rhs_ss] = ...
    convection_diffusion_1d( ...
    mesh,...
    U,...
    k,...
    bc,...
    @upwind_scheme_1d,...
    fid);

%% =========================================================================
% Explicit Euler
%% =========================================================================

logger(fid,'INFO','Running explicit');

[phi_exp,time_exp] = ...
    unsteady_convection_diffusion_1d( ...
    mesh,...
    U,...
    k,...
    bc,...
    @upwind_scheme_1d,...
    dt,...
    t_end,...
    'explicit_euler',...
    fid);

%% =========================================================================
% Implicit Euler
%% =========================================================================

logger(fid,'INFO','Running implicit');

[phi_imp,time_imp] = ...
    unsteady_convection_diffusion_1d( ...
    mesh,...
    U,...
    k,...
    bc,...
    @upwind_scheme_1d,...
    dt,...
    t_end,...
    'implicit_euler',...
    fid);

%% =========================================================================
% Validation 1
% Error relative to steady solution
%% =========================================================================

err_exp = zeros(length(time_exp),1);
err_imp = zeros(length(time_imp),1);

for n = 1:length(time_exp)

    err_exp(n) = ...
        sqrt(mean( ...
        (phi_exp(:,n)-phi_ss).^2));

    err_imp(n) = ...
        sqrt(mean( ...
        (phi_imp(:,n)-phi_ss).^2));

end

%% =========================================================================
% Validation 2
% Residual history
%% =========================================================================

res_exp = zeros(length(time_exp),1);
res_imp = zeros(length(time_imp),1);

for n = 1:length(time_exp)

    phi_e = phi_exp(2:end-1,n);
    phi_i = phi_imp(2:end-1,n);

    r_exp = ...
        A_ss.lower_diag.*[0;phi_e(1:end-1)] ...
      + A_ss.main_diag .*phi_e ...
      + A_ss.upper_diag.*[phi_e(2:end);0] ...
      - rhs_ss;

    r_imp = ...
        A_ss.lower_diag.*[0;phi_i(1:end-1)] ...
      + A_ss.main_diag .*phi_i ...
      + A_ss.upper_diag.*[phi_i(2:end);0] ...
      - rhs_ss;

    res_exp(n) = norm(r_exp,2);
    res_imp(n) = norm(r_imp,2);

end

%% =========================================================================
% Figure 1
% Error history
%% =========================================================================

f1 = figure('Color','w');

semilogy(time_exp,...
         err_exp,...
         'b',...
         'LineWidth',2);

hold on

semilogy(time_imp,...
         err_imp,...
         'r',...
         'LineWidth',2);

xlabel('t');
ylabel('L_2 Error');

title('Error Relative to Steady Solution');

legend('Explicit',...
       'Implicit');

grid on
box on

saveas(f1,...
    fullfile(figures_dir,...
    'error_history.png'));

%% =========================================================================
% Figure 2
% Residual history
%% =========================================================================

f2 = figure('Color','w');

semilogy(time_exp,...
         res_exp,...
         'b',...
         'LineWidth',2);

hold on

semilogy(time_imp,...
         res_imp,...
         'r',...
         'LineWidth',2);

xlabel('t');
ylabel('Residual Norm');

title('Residual Convergence');

legend('Explicit',...
       'Implicit');

grid on
box on

saveas(f2,...
    fullfile(figures_dir,...
    'residual_history.png'));

%% =========================================================================
% Figure 3
% Final solution
%% =========================================================================

f3 = figure('Color','w');

plot(x,...
     phi_ss,...
     'k--',...
     'LineWidth',3);

hold on

plot(x,...
     phi_exp(:,end),...
     'b',...
     'LineWidth',2);

plot(x,...
     phi_imp(:,end),...
     'r',...
     'LineWidth',2);

xlabel('x');
ylabel('\phi');

title('Final Solution');

legend('Steady',...
       'Explicit',...
       'Implicit');

grid on
box on

saveas(f3,...
    fullfile(figures_dir,...
    'final_solution.png'));

%% =========================================================================
% Figure 4
% Centre node history
%% =========================================================================

imid = round(length(x)/2);

f4 = figure('Color','w');

plot(time_exp,...
     phi_exp(imid,:),...
     'b',...
     'LineWidth',2);

hold on

plot(time_imp,...
     phi_imp(imid,:),...
     'r',...
     'LineWidth',2);

yline(phi_ss(imid),...
      'k--',...
      'LineWidth',2);

xlabel('t');
ylabel('\phi');

title('Centre Node History');

legend('Explicit',...
       'Implicit',...
       'Steady');

grid on
box on

saveas(f4,...
    fullfile(figures_dir,...
    'centre_node_history.png'));

%% =========================================================================
% Final statistics
%% =========================================================================

fprintf('\n');

fprintf('Explicit Final Error   = %.6e\n',err_exp(end));
fprintf('Implicit Final Error   = %.6e\n',err_imp(end));

fprintf('Explicit Final Residual = %.6e\n',res_exp(end));
fprintf('Implicit Final Residual = %.6e\n',res_imp(end));

fprintf('\n');

logger(fid,'INFO',...
    'Explicit Final Error = %.6e',err_exp(end));

logger(fid,'INFO',...
    'Implicit Final Error = %.6e',err_imp(end));

logger(fid,'INFO',...
    'Explicit Final Residual = %.6e',res_exp(end));

logger(fid,'INFO',...
    'Implicit Final Residual = %.6e',res_imp(end));

%% =========================================================================
% Save data
%% =========================================================================

logger(fid,'INFO','Simulation completed');

logger_close(fid,problem_name);
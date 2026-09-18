%% =========================================================================
% CFD - 1D Steady Convection Diffusion
% Test using B3 (Second-Order Upwind)
% =========================================================================

clear
clc
close all

%% Logging

problem_name = 'b3_1d_test';

fid = logger_init(problem_name);

logger(fid,'INFO','Starting simulation');

%% Parameters

L = 10.0;
U = 50.0;
k = 1e-3;

bc.W.type  = 'dirichlet';
bc.W.value = 0.0;

bc.E.type  = 'dirichlet';
bc.E.value = 20.0;

N = 1e3;

logger(fid,'INFO','L = %.3f',L);
logger(fid,'INFO','U = %.3f',U);
logger(fid,'INFO','k = %.6e',k);

logger(fid,'INFO','West BC : %s %.3f', ...
    bc.W.type,bc.W.value);

logger(fid,'INFO','East BC : %s %.3f', ...
    bc.E.type,bc.E.value);

logger(fid,'INFO','N = %d',N);

logger(fid,'INFO', ...
    'Scheme = B3 (2nd Order Upwind)');

%% Create mesh

logger(fid,'INFO','Creating uniform mesh');

mesh = create_uniform_mesh_1d(L,N);

%% Solve

logger(fid,'INFO', ...
    'Solving convection-diffusion equation');

[phi,A,rhs] = convection_diffusion_1d( ...
    mesh,...
    U,...
    k,...
    bc,...
    @b3_scheme_1d,...
    fid);

%% Plot

logger(fid,'INFO','Generating plot');

figure('Color','white');

plot(mesh.x,phi,'LineWidth',2);

xlabel('x');
ylabel('\phi');

title('1D Convection-Diffusion (B3 Scheme)');

grid on;
box on;

%% Create output folders

results_dir = fullfile('results',problem_name);

figures_dir = fullfile(results_dir,'figures');
data_dir    = fullfile(results_dir,'data');

if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

if ~exist(data_dir,'dir')
    mkdir(data_dir);
end

%% Save figure

saveas( ...
    gcf,...
    fullfile(figures_dir,'b3_solution.png'));

logger(fid,'INFO','Figure saved');

%% Save solution

save( ...
    fullfile(data_dir,'b3_solution.mat'),...
    'mesh',...
    'phi',...
    'A',...
    'rhs');

logger(fid,'INFO','Solution saved');

%% Finish

logger(fid,'INFO', 'Simulation completed successfully');

logger_close(fid,problem_name);
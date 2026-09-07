%% =========================================================================
% CFD - 1D Steady Convection Diffusion
% Test using Downwind discretization
% =========================================================================

clear
clc
close all

%% Logging

problem_name = 'downwind1d_test';

fid = logger_init(problem_name);

logger(fid,'INFO','Starting simulation');

%% Parameters

L = 10.0;
U = -5.0;
k = 1e-3;

phi_W = 0.0;
phi_E = 20.0;

N = 100;

logger(fid,'INFO','L = %.3f',L);
logger(fid,'INFO','U = %.3f',U);
logger(fid,'INFO','k = %.6e',k);
logger(fid,'INFO','phi_W = %.3f',phi_W);
logger(fid,'INFO','phi_E = %.3f',phi_E);
logger(fid,'INFO','N = %d',N);
logger(fid,'INFO','Scheme = Downwind');

%% Create mesh

logger(fid,'INFO','Creating uniform mesh');

mesh = create_uniform_mesh(L,N);

%% Solve

logger(fid,'INFO','Solving convection-diffusion equation');

[phi,A,rhs] = convection_diffusion_1d(mesh,U,k,phi_W,phi_E,@downwind_scheme,fid);

%% Plot

logger(fid,'INFO','Generating plot');

figure('Color','white');

plot(mesh.x,phi,'LineWidth',2);

xlabel('x');
ylabel('\phi');

title('1D Convection-Diffusion (Downwind Scheme)');

grid on;
box on;

%% Create output folders

results_dir = fullfile('results',problem_name);

figures_dir = fullfile(results_dir,'figures');
data_dir = fullfile(results_dir,'data');

if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

if ~exist(data_dir,'dir')
    mkdir(data_dir);
end

%% Save figure

saveas(gcf,fullfile(figures_dir,'downwind_solution.png'));

logger(fid,'INFO','Figure saved');

%% Save solution

save(fullfile(data_dir,'downwind_solution.mat'),'mesh','phi','A','rhs');

logger(fid,'INFO','Solution saved');

%% Finish

logger(fid,'INFO','Simulation completed successfully');

logger_close(fid,problem_name);
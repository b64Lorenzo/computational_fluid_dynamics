%% =========================================================================
% CFD - 1D Steady Convection Diffusion
% Test using Upwind discretization
% =========================================================================

clear
clc
close all

%% Logging

problem_name = 'upwind1d_test';

fid = logger_init(problem_name);

logger(fid,'INFO','Starting simulation');

%% Parameters

L = 10.0;
U = 50.0;
k = 1e-3;

bc.W.type = 'dirichlet';
bc.W.value = 0.0;

bc.E.type = 'dirichlet';
bc.E.value = 20.0;

N = 100;

logger(fid,'INFO','L = %.3f',L);
logger(fid,'INFO','U = %.3f',U);
logger(fid,'INFO','k = %.6e',k);
logger(fid,'INFO', 'WEST Bc type: ', bc.W.type);
logger(fid,'INFO','West BC  : %s %.3f',bc.W.type,bc.W.value);
logger(fid,'INFO', 'East Bc type: ', bc.E.type);
logger(fid,'INFO','East BC  : %s %.3f',bc.E.type,bc.E.value);
logger(fid,'INFO','N = %d',N);
logger(fid,'INFO','Scheme = Upwind');

%% Create Mesh

logger(fid,'INFO','Creating uniform mesh');

mesh = create_uniform_mesh_1d(L,N);

%% Solve

logger(fid,'INFO','Solving convection-diffusion equation');

[phi,A,rhs] = convection_diffusion_1d(mesh,U,k,bc,@upwind_scheme_1d,fid);
%% Plot

logger(fid,'INFO','Generating plot');

figure('Color','w');

plot(mesh.x,phi,'LineWidth',2);

xlabel('x');
ylabel('\phi');

title('1D Convection-Diffusion (Upwind Scheme)');

grid on;
box on;

%% Create Output Directories

results_dir = fullfile('results',problem_name);

figures_dir = fullfile(results_dir,'figures');
data_dir = fullfile(results_dir,'data');

if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

if ~exist(data_dir,'dir')
    mkdir(data_dir);
end

%% Save Figure

saveas(gcf,fullfile(figures_dir,'upwind_solution.png'));

logger(fid,'INFO','Figure saved');

%% Save Solution

save(fullfile(data_dir,'upwind_solution.mat'),'mesh','phi','A','rhs');

logger(fid,'INFO','Solution saved');

%% Finish

logger(fid,'INFO','Simulation completed successfully');

logger_close(fid,problem_name);

%% =========================================================================
% CFD - 1D Steady Convection Diffusion
% Test using Method B on a Non-Uniform Grid
% =========================================================================

clear
clc
close all

%% Logging

problem_name = 'methodB_1d_test';

fid = logger_init(problem_name);

logger(fid,'INFO','Starting simulation');

%% Parameters

L = 1.0;
U = 1.0;
k = 2e-3;

bc.W.type  = 'dirichlet';
bc.W.value = 0.0;

bc.E.type  = 'dirichlet';
bc.E.value = 1.0;

N = 1e4;

logger(fid,'INFO','L = %.3f',L);
logger(fid,'INFO','U = %.3f',U);
logger(fid,'INFO','k = %.6e',k);
logger(fid,'INFO','N = %d',N);

logger(fid,'INFO', 'Scheme = Method A');

%% Create mesh

logger(fid,'INFO', ...
    'Creating non-uniform mesh');

mesh = nonuniform_121(L,k);

%% Solve

logger(fid,'INFO', ...
    'Solving convection-diffusion equation');

[phi,A,rhs] = convection_diffusion_1d( ...
    mesh,...
    U,...
    k,...
    bc,...
    @convection_methodB_1d,...
    fid);

%% Plot solution

logger(fid,'INFO', ...
    'Generating plots');

figure('Color','white');

plot(mesh.x,phi,'LineWidth',2);

xlabel('x');
ylabel('\phi');

title('1D Convection-Diffusion (Method B)');

grid on;
box on;

%% Plot grid spacing

figure('Color','white');

plot(mesh.x(1:end-1),diff(mesh.x), ...
    'LineWidth',2);

xlabel('x');
ylabel('\Delta x');

title('Non-uniform Grid Spacing');

grid on;
box on;

%% Output folders

results_dir = fullfile( ...
    'results',problem_name);

figures_dir = fullfile( ...
    results_dir,'figures');

data_dir = fullfile( ...
    results_dir,'data');

if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

if ~exist(data_dir,'dir')
    mkdir(data_dir);
end

%% Save figures

figure(1);

saveas(gcf,...
    fullfile(figures_dir,...
    'methodB_solution.png'));

figure(2);

saveas(gcf,...
    fullfile(figures_dir,...
    'methodB_grid.png'));

%% Save data

save( ...
    fullfile(data_dir,...
    'methodB_solution.mat'), ...
    'mesh',...
    'phi',...
    'A',...
    'rhs');

logger(fid,'INFO', ...
    'Solution saved');

%% Finish

logger(fid,'INFO', ...
    'Simulation completed successfully');

logger_close(fid,problem_name);
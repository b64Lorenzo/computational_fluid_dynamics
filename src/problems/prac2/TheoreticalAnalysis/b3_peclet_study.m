%% =========================================================================
% CFD - B3 Wiggle-Free Verification
%
% Verify numerically that the B3 scheme remains wiggle-free for all mesh
% Peclet numbers by varying the grid size.
%
% Parameters:
%   U = 5
%   k = 0.01
%   L = 1
%   phi(0) = 0
%   phi(L) = 1
%
% The mesh Peclet number is
%
%   P = U*h/k = U*L/(k*N)
%
%========================================================================

clear
clc
close all

%% Logging

problem_name = 'b3_peclet_study';

fid = logger_init(problem_name);

logger(fid,'INFO','Starting simulation');

%% Parameters

L = 1.0;
U = 5.0;
k = 0.01;

bc.W.type  = 'dirichlet';
bc.W.value = 0.0;

bc.E.type  = 'dirichlet';
bc.E.value = 1.0;

logger(fid,'INFO','L = %.3f',L);
logger(fid,'INFO','U = %.3f',U);
logger(fid,'INFO','k = %.6e',k);

%% Grid cases

N_list = [5 10 20 50 100];

logger(fid,'INFO','Running %d grid cases',length(N_list));

%% Plot

figure('Color','white');
hold on

legend_entries = cell(length(N_list),1);

%% Loop over grids

for icase = 1:length(N_list)

    N = N_list(icase);

    logger(fid,'INFO', '------------------------------------------------');

    logger(fid,'INFO','Case %d/%d',icase,length(N_list));

    logger(fid,'INFO','N = %d',N);

    %% Mesh

    mesh = create_uniform_mesh_1d(L,N);

    h = L/N;

    P = U*h/k;

    logger(fid,'INFO', 'h = %.6e',h);

    logger(fid,'INFO', 'Mesh Peclet number = %.6f',P);

    %% Solve

    [phi,A,rhs] = convection_diffusion_1d( ...
        mesh,...
        U,...
        k,...
        bc,...
        @b3_scheme_1d,...
        fid);

    %% Plot

    plot(mesh.x,phi,'LineWidth',2);

    legend_entries{icase} = sprintf('N=%d, P=%.1f',N,P);

    %% Save case

    results_dir = fullfile('results',problem_name);

    data_dir = fullfile( results_dir,'data');

    if ~exist(data_dir,'dir')
        mkdir(data_dir);
    end

    save( ...
        fullfile(data_dir,...
        sprintf('case_N%d.mat',N)),'mesh','phi','A','rhs','P');

end

%% Finish figure

xlabel('x');
ylabel('\phi');

title('B3 Scheme for Different Mesh Peclet Numbers');

legend(legend_entries,'Location','best');

grid on;
box on;

%% Create folders

results_dir = fullfile( 'results',problem_name);

figures_dir = fullfile( results_dir,'figures');

if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

%% Save figure

saveas(gcf,fullfile(figures_dir,'b3_peclet_study.png'));

logger(fid,'INFO', 'Figure saved');

%% Finish

logger(fid,'INFO', 'Study completed successfully');

logger_close(fid,problem_name);
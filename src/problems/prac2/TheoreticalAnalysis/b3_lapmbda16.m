%% =========================================================================
% Compare B3 and lambda = 1/6
% One figure
% =========================================================================

clear
clc
close all

%% Logging

problem_name = 'b3_vs_lambda16';

fid = logger_init(problem_name);

%% Parameters

L = 1.0;
U = 5.0;
k = 1e-6;

bc.W.type  = 'dirichlet';
bc.W.value = 0.0;

bc.E.type  = 'dirichlet';
bc.E.value = 1.0;

%% Grids

N_list = [25 50 100];

%% Figure

figure('Color','white');
hold on

colors = lines(length(N_list));

legend_entries = {};

%% Loop

for icase = 1:length(N_list)

    N = N_list(icase);

    mesh = create_uniform_mesh_1d(L,N);

    h = L/N;

    P = U*h/k;

    logger(fid,'INFO','N=%d  P=%e',N,P);

    %% B3

    [phi_b3,~,~] = convection_diffusion_1d( ...
        mesh,...
        U,...
        k,...
        bc,...
        @b3_scheme_1d,...
        fid);

    %% lambda = 1/6

    lambda16 = ...
        @(mesh,U,k,i) ...
        lambda_scheme_1d( ...
        mesh,U,k,i,1/6);

    [phi_l16,~,~] = convection_diffusion_1d( ...
        mesh,...
        U,...
        k,...
        bc,...
        lambda16,...
        fid);

    %% Plot

    plot(mesh.x,...
        phi_b3,...
        '-',...
        'LineWidth',2,...
        'Color',colors(icase,:));

    plot(mesh.x,...
        phi_l16,...
        '--',...
        'LineWidth',2,...
        'Color',colors(icase,:));

    legend_entries{end+1} = sprintf('B3, N=%d, P=%.1e',N,P);

    legend_entries{end+1} =sprintf('\\lambda=1/6, N=%d, P=%.1e',N,P);

end

xlabel('x');
ylabel('\phi');

title('B3 versus \lambda = 1/6');

legend( legend_entries,'Location','best');

grid on
box on

%% Save

results_dir = fullfile('results',problem_name);

figures_dir = fullfile(results_dir,'figures');

if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

saveas( ...
    gcf,...
    fullfile(figures_dir,...
    'b3_vs_lambda16.png'));

logger_close(fid,problem_name);
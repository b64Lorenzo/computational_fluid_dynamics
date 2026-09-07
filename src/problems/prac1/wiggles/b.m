%% =========================================================================
% CFD - Exercise 2(b)
%
% Verify:
%
%       delta ~ h
%
% for the upwind scheme.
%
% =========================================================================

clear
clc
close all

%% Logging

problem_name = 'prac1_upwind_boundary_layer_ex2b';

fid = logger_init(problem_name);

try

    logger(fid,'INFO','Starting simulation');

    %% Parameters

    L = 10.0;
    U = 50.0;
    k = 1e-3;

    bc.W.type  = 'dirichlet';
    bc.W.value = 0.0;

    bc.E.type  = 'dirichlet';
    bc.E.value = 20.0;

    %% Meshes

    N_values = [ ...
        100 ...
        250 ...
        500 ...
        1000 ...
        2500 ...
        5000 ...
        10000 ...
        25000 ...
        50000 ...
        100000 ...
        250000 ...
        500000 ...
        1000000];

    delta_num = zeros(size(N_values));

    %% Loop

    for idx = 1:length(N_values)

        N = N_values(idx);

        logger(fid,'INFO',...
            'Running case N = %d',N);

        mesh = create_uniform_mesh_1d(L,N);

        [phi,~,~] = convection_diffusion_1d( ...
            mesh,...
            U,...
            k,...
            bc,...
            @upwind_scheme_1d,...
            fid);

        x = mesh.x;

        %% Boundary-layer thickness

        phi_mid = 0.5*( ...
            bc.W.value + bc.E.value);

        idx_cross = find(phi >= phi_mid,1,'first');

        if isempty(idx_cross)

            delta_num(idx) = NaN;

        elseif idx_cross == 1

            delta_num(idx) = NaN;

        else

            x1   = x(idx_cross-1);
            x2   = x(idx_cross);

            phi1 = phi(idx_cross-1);
            phi2 = phi(idx_cross);

            % Linear interpolation

            x50 = x1 + ...
                 (phi_mid-phi1) * ...
                 (x2-x1)/(phi2-phi1);

            delta_num(idx) = L - x50;

        end

        logger(fid,'INFO',...
            'N=%d | delta=%.6e',...
            N,...
            delta_num(idx));

    end

    %% Remove NaN values

    valid_idx = ~isnan(delta_num);

    N_plot     = N_values(valid_idx);
    delta_plot = delta_num(valid_idx);

    h_plot = L ./ N_plot;

    %% Fit delta = C*h^p

    p = polyfit( ...
        log10(h_plot), ...
        log10(delta_plot), ...
        1);

    observed_slope = p(1);

    logger(fid,'INFO',...
        'Observed slope = %.3f',...
        observed_slope);

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

    %% =====================================================================
    % Figure : delta vs h
    % =====================================================================

    figure('Color','white');

    loglog(h_plot,...
           delta_plot,...
           'bo-',...
           'LineWidth',2,...
           'MarkerSize',8,...
           'MarkerFaceColor','b');

    hold on

    C = 10^(p(2));

    loglog(h_plot,...
           C*h_plot.^observed_slope,...
           'r--',...
           'LineWidth',2);

    xlabel('Grid Size, h');
    ylabel('\delta');

    title({ ...
        'Boundary Layer Thickness vs Grid Size', ...
        sprintf('Observed slope = %.3f', ...
        observed_slope)});

    legend( ...
        'Numerical',...
        sprintf('\\delta \\propto h^{%.2f}',...
        observed_slope),...
        'Location','northwest');

    grid on
    box on

    saveas(gcf,...
        fullfile(figures_dir,...
        'boundary_layer_vs_h_loglog.png'));

    %% Save

    save( ...
        fullfile(data_dir,...
        'boundary_layer_scaling.mat'), ...
        'N_values',...
        'delta_num',...
        'observed_slope');

    logger(fid,'INFO',...
        'Simulation completed successfully');

catch ME

    logger(fid,'ERROR',...
        'Fatal error occurred.');

    logger_exception(fid,ME);

end

logger_close(fid,problem_name);
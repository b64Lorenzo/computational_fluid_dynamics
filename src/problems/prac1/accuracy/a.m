%% =========================================================================
% CFD - Exercise 3
%
% Comparison of Central and Upwind Differencing
%
% Pe = 500
%
% Order from:
% 1) Analytical solution
% 2) Consecutive meshes (i and i+1)
%
% =========================================================================

clear
clc
close all

%% Logging

problem_name = 'prac1_wiggles_ex3';

fid = logger_init(problem_name);

try

    logger(fid,'INFO','Starting simulation');

    %% Parameters

    L = 10.0;
    U = 50.0;
    k = 1.0;

    bc.W.type  = 'dirichlet';
    bc.W.value = 0.0;

    bc.E.type  = 'dirichlet';
    bc.E.value = 20.0;

    Pe = U*L/k;

    logger(fid,'INFO','Global Peclet number = %.2f',Pe);

    %% Meshes

    N_values = [ ...
        5000 ...
        10000 ...
        25000 ...
        50000 ...
        100000 ...
        200000 ...
        250000 ...
        500000 ...
        1000000];

    %% Output folders

    results_dir = fullfile('results',problem_name);

    figures_dir = fullfile(results_dir,'figures');
    data_dir    = fullfile(results_dir,'data');

    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    if ~exist(data_dir,'dir')
        mkdir(data_dir);
    end

    %% =====================================================================
    % Compute solutions
    %% =====================================================================

    central_solutions = cell(length(N_values),1);
    upwind_solutions  = cell(length(N_values),1);
    meshes            = cell(length(N_values),1);

    errors_exact_central = zeros(size(N_values));
    errors_exact_upwind  = zeros(size(N_values));

    for idx = 1:length(N_values)

        N = N_values(idx);

        logger(fid,'INFO','Running N = %d',N);

        mesh = create_uniform_mesh_1d(L,N);

        meshes{idx} = mesh;

        %% Central

        [phi_central,~,~] = convection_diffusion_1d( ...
            mesh,...
            U,...
            k,...
            bc,...
            @central_scheme_1d,...
            []);

        %% Upwind

        [phi_upwind,~,~] = convection_diffusion_1d( ...
            mesh,...
            U,...
            k,...
            bc,...
            @upwind_scheme_1d,...
            []);

        central_solutions{idx} = phi_central;
        upwind_solutions{idx}  = phi_upwind;

        x = mesh.x;

        %% Stable analytical solution

        phi_exact = bc.W.value + ...
            (bc.E.value-bc.W.value).* ...
            (exp(-Pe) - exp(-Pe*(1-x/L))) ./ ...
            (exp(-Pe)-1);

        %% Error vs exact

        errors_exact_central(idx) = ...
            norm(phi_central-phi_exact,2) ...
            / sqrt(length(phi_exact));

        errors_exact_upwind(idx) = ...
            norm(phi_upwind-phi_exact,2) ...
            / sqrt(length(phi_exact));

    end

    %% =====================================================================
    % Orders from analytical solution
    %% =====================================================================

    order_exact_central = zeros(length(N_values)-1,1);
    order_exact_upwind  = zeros(length(N_values)-1,1);

    for i = 1:length(order_exact_central)

        h1 = L/N_values(i);
        h2 = L/N_values(i+1);

        order_exact_central(i) = ...
            log(errors_exact_central(i+1) ...
            / errors_exact_central(i)) ...
            / log(h2/h1);

        order_exact_upwind(i) = ...
            log(errors_exact_upwind(i+1) ...
            / errors_exact_upwind(i)) ...
            / log(h2/h1);

    end

    avg_order_exact_central = mean(order_exact_central);
    avg_order_exact_upwind  = mean(order_exact_upwind);

    %% =====================================================================
    % Consecutive-mesh errors
    %% =====================================================================

    errors_self_central = zeros(length(N_values)-1,1);
    errors_self_upwind  = zeros(length(N_values)-1,1);

    for i = 1:length(N_values)-1

        x_coarse = meshes{i}.x;

        %% Central

        phi_fine_interp = interp1( ...
            meshes{i+1}.x,...
            central_solutions{i+1},...
            x_coarse,...
            'linear');

        errors_self_central(i) = ...
            norm(central_solutions{i}-phi_fine_interp,2) ...
            / sqrt(length(x_coarse));

        %% Upwind

        phi_fine_interp = interp1( ...
            meshes{i+1}.x,...
            upwind_solutions{i+1},...
            x_coarse,...
            'linear');

        errors_self_upwind(i) = ...
            norm(upwind_solutions{i}-phi_fine_interp,2) ...
            / sqrt(length(x_coarse));

    end

    %% =====================================================================
    % Orders from consecutive meshes
    %% =====================================================================

    order_self_central = zeros(length(errors_self_central)-1,1);
    order_self_upwind  = zeros(length(errors_self_upwind)-1,1);

    for i = 1:length(order_self_central)

        h1 = L/N_values(i);
        h2 = L/N_values(i+1);

        order_self_central(i) = ...
            log(errors_self_central(i+1) ...
            / errors_self_central(i)) ...
            / log(h2/h1);

        order_self_upwind(i) = ...
            log(errors_self_upwind(i+1) ...
            / errors_self_upwind(i)) ...
            / log(h2/h1);

    end

    avg_order_self_central = mean(order_self_central);
    avg_order_self_upwind  = mean(order_self_upwind);

    %% =====================================================================
    % Figure 1 : Analytical solution errors
    %% =====================================================================

    figure('Color','white');

    loglog(N_values,...
           errors_exact_central,...
           'bo-',...
           'LineWidth',2,...
           'MarkerFaceColor','b');

    hold on

    loglog(N_values,...
           errors_exact_upwind,...
           'rs-',...
           'LineWidth',2,...
           'MarkerFaceColor','r');

    xlabel('Number of Cells, N');
    ylabel('L_2 Error');

    title({ ...
        'Error Relative to Analytical Solution', ...
        sprintf(['Central = %.2f | Upwind = %.2f'], ...
        avg_order_exact_central,...
        avg_order_exact_upwind)});

    legend('Central','Upwind',...
        'Location','southwest');

    grid on
    box on

    saveas(gcf,...
        fullfile(figures_dir,...
        'error_vs_analytical.png'));

    %% =====================================================================
    % Figure 2 : Consecutive-mesh convergence
    %% =====================================================================

    figure('Color','white');

    loglog(N_values(1:end-1),...
           errors_self_central,...
           'bo-',...
           'LineWidth',2,...
           'MarkerFaceColor','b');

    hold on

    loglog(N_values(1:end-1),...
           errors_self_upwind,...
           'rs-',...
           'LineWidth',2,...
           'MarkerFaceColor','r');

    xlabel('Number of Cells, N');
    ylabel('L_2 Difference');

    title({ ...
        'Convergence from Consecutive Meshes', ...
        sprintf(['Central = %.2f | Upwind = %.2f'], ...
        avg_order_self_central,...
        avg_order_self_upwind)});

    legend('Central','Upwind',...
        'Location','southwest');

    grid on
    box on

    saveas(gcf,...
        fullfile(figures_dir,...
        'self_convergence.png'));

    %% Logging

    logger(fid,...
        'INFO',...
        'Central order vs exact = %.3f',...
        avg_order_exact_central);

    logger(fid,...
        'INFO',...
        'Upwind order vs exact = %.3f',...
        avg_order_exact_upwind);

    logger(fid,...
        'INFO',...
        'Central order from consecutive meshes = %.3f',...
        avg_order_self_central);

    logger(fid,...
        'INFO',...
        'Upwind order from consecutive meshes = %.3f',...
        avg_order_self_upwind);

    %% Save

    save( ...
        fullfile(data_dir,'exercise3_convergence.mat'), ...
        'N_values',...
        'errors_exact_central',...
        'errors_exact_upwind',...
        'errors_self_central',...
        'errors_self_upwind',...
        'order_exact_central',...
        'order_exact_upwind',...
        'order_self_central',...
        'order_self_upwind');

    logger(fid,'INFO','Data saved');
    logger(fid,'INFO','Simulation completed successfully');

catch ME

    logger(fid,'ERROR','Fatal error occurred.');
    logger_exception(fid,ME);

end

logger_close(fid,problem_name);
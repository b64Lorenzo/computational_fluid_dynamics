%% =========================================================================
% CFD - Exercise 2(c)
% Wiggles Study using Central Differencing
%
% Investigate when central differencing becomes wiggle-free.
% =========================================================================

clear
clc
close all

%% Logging

problem_name = 'prac1_wiggles_ex2c';

fid = logger_init(problem_name);

try

    logger(fid,'INFO','Starting simulation');

    %% Parameters

    L = 10.0;
    U = 50.0;
    k = 1e-3;

    bc.W.type  = 'dirichlet';
    bc.W.value = 0.0;

    bc.E.type  = 'neumann';
    bc.E.value = 20;

    %% Theoretical limit

    N_crit = ceil(U*L/(2*k));

    logger(fid,'INFO','Critical cell Peclet criterion: Pe_h <= 2');
    logger(fid,'INFO','Theoretical critical N = %d',N_crit);

    %% Mesh sizes to investigate

    N_values = [5000, 10000, 25000 50000 100000 200000 225000 250000 500000 1000000];

    Pe_values = zeros(size(N_values));

    %% Figure 1 : Solutions

    figure('Color','white');

    hold on

    colors = lines(length(N_values));

    successful_cases = [];

    solutions = cell(length(N_values),1);
    meshes = cell(length(N_values),1);

    for idx = 1:length(N_values)

        N = N_values(idx);

        try

            logger(fid,'INFO','Running case N = %d',N);

            mesh = create_uniform_mesh_1d(L,N);

            [phi,~,~] = convection_diffusion_1d( ...
                mesh,U,k,bc,@central_scheme_1d,fid);
            
            solutions{idx} = phi;
            meshes{idx} = mesh;

            dx = L/N;

            Pe_h = U*dx/k;

            Pe_values(idx) = Pe_h;

            successful_cases(end+1) = N;

            logger(fid,'INFO',...
                'N = %d | dx = %.6e | Pe_h = %.3f',...
                N,dx,Pe_h);

            if Pe_h <= 2

                plot(mesh.x,phi,...
                    'Color',colors(idx,:),...
                    'LineWidth',3,...
                    'DisplayName',...
                    sprintf('N=%d (Pe_h=%.2f)',N,Pe_h));

            else

                plot(mesh.x,phi,...
                    '--',...
                    'Color',colors(idx,:),...
                    'LineWidth',1.5,...
                    'DisplayName',...
                    sprintf('N=%d (Pe_h=%.1f)',N,Pe_h));

            end

        catch ME_case

            logger(fid,'ERROR','Case N=%d failed.',N);

            logger_exception(fid,ME_case);

            continue

        end

    end

    xlabel('x');
    ylabel('\phi');

    title('Exercise 2(a): Central Differencing Solutions');

    legend('Location','eastoutside');

    grid on;
    box on;

    xlim([9.75 10.0]);

    %% Annotation

    txt = sprintf( ...
        'Theoretical limit: Pe_h <= 2\nCritical N = %d', ...
        N_crit);

    annotation('textbox',...
        [0.15 0.75 0.25 0.10],...
        'String',txt,...
        'FitBoxToText','on',...
        'BackgroundColor','white');

    %% Save Figure 1

    results_dir = fullfile('results',problem_name);

    figures_dir = fullfile(results_dir,'figures');
    data_dir    = fullfile(results_dir,'data');

    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    if ~exist(data_dir,'dir')
        mkdir(data_dir);
    end

    saveas(gcf,...
        fullfile(figures_dir,'wiggles_study_central.png'));

    logger(fid,'INFO','Solution figure saved');

    %% =====================================================================
        %% =====================================================================
    % Figure 2 : Convergence Study Relative to Reference Solution
    % =====================================================================
    
    reference_idx = find(N_values == 1000000,1);
    
    phi_ref = solutions{reference_idx};
    x_ref   = meshes{reference_idx}.x;
    
    errors = NaN(size(N_values));
    
    for idx = 1:length(N_values)
    
        if idx == reference_idx
            continue
        end
    
        phi_current = solutions{idx};
        x_current   = meshes{idx}.x;
    
        phi_ref_interp = interp1( ...
            x_ref,...
            phi_ref,...
            x_current,...
            'linear');
    
        errors(idx) = norm(phi_current - phi_ref_interp,2) ...
                    / sqrt(length(phi_current));
    
    end
    
    %% ---------------------------------------------------------
    % Local convergence slopes
    % ---------------------------------------------------------
    
    valid_idx = find(~isnan(errors));
    
    N_fit   = N_values(valid_idx);
    err_fit = errors(valid_idx);
    
    local_slope = zeros(length(err_fit)-1,1);
    
    for i = 1:length(local_slope)
    
        local_slope(i) = ...
            (log(err_fit(i+1)) - log(err_fit(i))) / ...
            (log(N_fit(i+1))   - log(N_fit(i)));
    
    end
    
    %% ---------------------------------------------------------
    % Average slope before / after Pe_h = 2
    % ---------------------------------------------------------
    
    slope_midpoints = sqrt(N_fit(1:end-1).*N_fit(2:end));
    
    before_idx = slope_midpoints < N_crit;
    after_idx  = slope_midpoints >= N_crit;
    
    avg_slope_before = mean(local_slope(before_idx));
    avg_slope_after  = mean(local_slope(after_idx));
    
    logger(fid,'INFO', ...
        'Average slope before Pe_h=2 : %.3f', ...
        avg_slope_before);
    
    logger(fid,'INFO', ...
        'Average slope after  Pe_h=2 : %.3f', ...
        avg_slope_after);
    
    %% ---------------------------------------------------------
    % Find first interval with slope approximately -2
    % ---------------------------------------------------------
    
    tol = 0.2;
    
    idx_slope2 = find(abs(local_slope + 2) < tol,1,'first');
    
    highlight_idx = [];
    
    if ~isempty(idx_slope2)
    
        highlight_idx = idx_slope2 + 1;
    
        logger(fid,...
            'INFO',...
            'First slope near -2 at N = %d (slope = %.3f)',...
            N_fit(highlight_idx),...
            local_slope(idx_slope2));
    
    end
    
    %% ---------------------------------------------------------
    % Plot
    % ---------------------------------------------------------
    
    figure('Color','white');
    
    hold on
    
    loglog(N_fit,...
           err_fit,...
           'bo-',...
           'LineWidth',2,...
           'MarkerSize',8,...
           'MarkerFaceColor','b');
    
    %% Vertical line at Pe_h = 2
    
    xline(N_crit,...
        'r--',...
        sprintf('Pe_h = 2  (N = %d)',N_crit),...
        'LineWidth',2);
    
    %% Highlight first slope ≈ -2
    
    if ~isempty(highlight_idx)
    
        loglog( ...
            N_fit(highlight_idx), ...
            err_fit(highlight_idx), ...
            'kp', ...
            'MarkerSize',18, ...
            'MarkerFaceColor','y');
    
        text( ...
            1.05*N_fit(highlight_idx), ...
            err_fit(highlight_idx), ...
            sprintf('Slope = %.2f',local_slope(idx_slope2)), ...
            'FontWeight','bold');
    
    end
    
    %% Show local slope values
    
    for i = 1:length(local_slope)
    
        x_mid = sqrt(N_fit(i)*N_fit(i+1));
        y_mid = sqrt(err_fit(i)*err_fit(i+1));
    
        text( ...
            x_mid,...
            y_mid,...
            sprintf('%.2f',local_slope(i)),...
            'FontSize',8);
    
    end
    
    xlabel('Number of Cells, N');
    ylabel('L_2 Error');
    
    title({ ...
        'Convergence Relative to Reference Solution (N = 10^6)', ...
        sprintf(['Average slope before Pe_h=2: %.3f    |    ' ...
                  'Average slope after Pe_h=2: %.3f'], ...
                  avg_slope_before,...
                  avg_slope_after)});
    
    grid on
    box on
    
    annotation('textbox',...
        [0.14 0.73 0.34 0.12],...
        'String',sprintf([ ...
        'Pe_h = UL/(kN) \\leq 2\n' ...
        'Critical mesh: N = %d'], ...
        N_crit),...
        'FitBoxToText','on',...
        'BackgroundColor','white');
    
    saveas(gcf,...
        fullfile(figures_dir,'error_vs_reference_loglog.png'));
    
    logger(fid,'INFO',...
        'Error convergence figure saved');


    %% Save Data

    save( ...
        fullfile(data_dir,'wiggles_study_central.mat'), ...
        'N_values',...
        'Pe_values',...
        'successful_cases',...
        'L',...
        'U',...
        'k',...
        'N_crit');

    logger(fid,'INFO','Data saved');

    logger(fid,'INFO',...
        'Successful cases = %d/%d',...
        numel(successful_cases),...
        numel(N_values));

    logger(fid,'INFO','Simulation completed successfully');

catch ME

    logger(fid,'ERROR','Fatal error occurred.');

    logger_exception(fid,ME);

end

logger_close(fid,problem_name);
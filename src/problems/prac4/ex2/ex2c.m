%% =========================================================================
% EXERCISE 2(c)
%
% CENTRAL DIFFERENCING + EXPLICIT EULER
%
% Stability Investigation
%
% Initial condition:
%
%     phi(x,0) = 2 + sin(pi*x/L)
%
% Periodic boundary conditions:
%
%     phi(0) = phi(L)
%
% Investigate:
%
%     eta < d < 1
%
% versus
%
%     eta^2 < d < 1
%
%% =========================================================================

clear
clc
close all

%% =========================================================================
% Logging
%% =========================================================================

problem_name = 'exercise2c_central_explicit';

fid = logger_init(problem_name);

logger(fid,'INFO','Starting Exercise 2(c)');

%% =========================================================================
% Directories
%% =========================================================================

results_dir = fullfile('results',problem_name);

figures_dir = fullfile(results_dir,'figures');
data_dir    = fullfile(results_dir,'data');

if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

if ~exist(data_dir,'dir')
    mkdir(data_dir);
end

%% =========================================================================
% Parameters
%% =========================================================================

L = 1.0;

U = 1.0;

N = 101;

t_end = 0.05;

%% =========================================================================
% Uniform Mesh
%% =========================================================================

mesh = create_uniform_mesh_1d(L,N);

x = mesh.x(:);

dx = x(2)-x(1);

%% =========================================================================
% Boundary Conditions
%% =========================================================================

bc.W.type = 'periodic';
bc.E.type = 'periodic';

%% =========================================================================
% Initial Condition
%% =========================================================================

bc.IC = 2 + sin(pi*x/L);

%% =========================================================================
% Steady Solution
%% =========================================================================

phi_ss = 2*ones(size(x));

logger(fid,'INFO',...
    'Steady solution = phi = 2');

%% =========================================================================
% Figure: Initial Condition
%% =========================================================================

f0 = figure('Color','w');

plot(x,...
     bc.IC,...
     'LineWidth',2)

hold on

plot(x,...
     phi_ss,...
     'k--',...
     'LineWidth',2)

xlabel('x')
ylabel('\phi')

title('Initial Condition and Steady Solution')

legend('Initial condition',...
       'Steady solution',...
       'Location','best')

grid on
box on

saveas(f0,...
    fullfile(figures_dir,...
    'initial_condition.png'));

%% =========================================================================
% Parameter Space
%% =========================================================================

eta_values = [ ...
    0.10 ...
    0.20 ...
    0.30 ...
    0.40 ...
    0.50 ...
    0.60 ...
    0.80 ];

d_values = [ ...
    0.05 ...
    0.10 ...
    0.20 ...
    0.30 ...
    0.40 ...
    0.50 ...
    0.60 ...
    0.80 ...
    0.95 ];

nEta = length(eta_values);

nD = length(d_values);

%% =========================================================================
% Storage
%% =========================================================================

growth_map      = zeros(nEta,nD);

stable_map      = zeros(nEta,nD);

pos_coeff_map   = zeros(nEta,nD);

fourier_map     = zeros(nEta,nD);

%% =========================================================================
% Theoretical Predictions
%% =========================================================================

for iEta = 1:nEta

    eta = eta_values(iEta);

    for iD = 1:nD

        d = d_values(iD);

        %
        % Positive coefficient:
        %
        % eta < d < 1
        %

        if (eta < d) && (d < 1)

            pos_coeff_map(iEta,iD) = 1;

        end

        %
        % Fourier:
        %
        % eta^2 < d < 1
        %

        if (eta^2 < d) && (d < 1)

            fourier_map(iEta,iD) = 1;

        end

    end

end

%% =========================================================================
% Numerical Investigation
%% =========================================================================

ref_norm = norm(bc.IC-2,2);

for iEta = 1:nEta

    eta = eta_values(iEta);

    dt = eta*dx/U;

    logger(fid,'INFO',...
        'Testing eta = %.3f',...
        eta);

    for iD = 1:nD

        d = d_values(iD);

        %
        % d = k*dt/dx^2
        %

        k_eff = d*dx^2/dt;

        try

            [phi,time] = ...
                unsteady_convection_diffusion_1d( ...
                mesh,...
                U,...
                k_eff,...
                bc,...
                @central_scheme_1d,...
                dt,...
                t_end,...
                'explicit_euler',...
                fid);

            growth = ...
                norm(phi(:,end)-2,2) ...
               /norm(phi(:,1)-2,2);

            growth_map(iEta,iD) = growth;

            %
            % Stability criterion
            %

            if isfinite(growth) && growth < 2

                stable_map(iEta,iD) = 1;

            else

                stable_map(iEta,iD) = 0;

            end

        catch

            growth_map(iEta,iD) = Inf;

            stable_map(iEta,iD) = 0;

        end

    end

end

%% =========================================================================
% Figure 1
% Numerical Stability Map
%% =========================================================================

f1 = figure('Color','w');

imagesc(d_values,...
        eta_values,...
        stable_map)

axis xy

xlabel('d')
ylabel('\eta')

title('Numerical Stability Map')

colorbar

saveas(f1,...
    fullfile(figures_dir,...
    'numerical_stability_map.png'));

%% =========================================================================
% Figure 2
% Positive Coefficient Prediction
%% =========================================================================

f2 = figure('Color','w');

imagesc(d_values,...
        eta_values,...
        pos_coeff_map)

axis xy

xlabel('d')
ylabel('\eta')

title('Positive Coefficient Prediction')

colorbar

saveas(f2,...
    fullfile(figures_dir,...
    'positive_coefficient_prediction.png'));

%% =========================================================================
% Figure 3
% Fourier Prediction
%% =========================================================================

f3 = figure('Color','w');

imagesc(d_values,...
        eta_values,...
        fourier_map)

axis xy

xlabel('d')
ylabel('\eta')

title('Fourier Prediction')

colorbar

saveas(f3,...
    fullfile(figures_dir,...
    'fourier_prediction.png'));

%% =========================================================================
% Figure 4
% Amplification Map
%% =========================================================================

f4 = figure('Color','w');

imagesc(d_values,...
        eta_values,...
        log10(growth_map))

axis xy

xlabel('d')
ylabel('\eta')

title('log_{10}(Amplification Factor)')

colorbar

saveas(f4,...
    fullfile(figures_dir,...
    'amplification_map.png'));

%% =========================================================================
% Figure 5
% Comparison
%% =========================================================================

f5 = figure('Color','w');

subplot(1,3,1)

imagesc(d_values,...
        eta_values,...
        stable_map)

axis xy

title('Numerical')

xlabel('d')
ylabel('\eta')

subplot(1,3,2)

imagesc(d_values,...
        eta_values,...
        pos_coeff_map)

axis xy

title('\eta < d < 1')

xlabel('d')

subplot(1,3,3)

imagesc(d_values,...
        eta_values,...
        fourier_map)

axis xy

title('\eta^2 < d < 1')

xlabel('d')

saveas(f5,...
    fullfile(figures_dir,...
    'comparison.png'));

%% =========================================================================
% Save Results
%% =========================================================================

save(fullfile(data_dir,...
    'exercise2c_results.mat'),...
    'eta_values',...
    'd_values',...
    'stable_map',...
    'growth_map',...
    'pos_coeff_map',...
    'fourier_map');

%% =========================================================================
% Summary
%% =========================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('EXERCISE 2(c)\n');
fprintf('====================================================\n');
fprintf('\n');
fprintf('Steady solution:\n');
fprintf('phi_ss = 2\n');
fprintf('\n');
fprintf('Compare:\n');
fprintf('  Numerical Stability Map\n');
fprintf('  eta < d < 1\n');
fprintf('  eta^2 < d < 1\n');
fprintf('\n');
fprintf('The prediction that best matches\n');
fprintf('the numerical stability map is the\n');
fprintf('correct stability requirement.\n');
fprintf('====================================================\n');

logger(fid,'INFO','Exercise completed');

logger_close(fid,problem_name);
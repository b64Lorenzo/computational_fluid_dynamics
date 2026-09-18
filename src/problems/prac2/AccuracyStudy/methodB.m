%% =========================================================================
% CFD - Method B
%
% Accuracy and Eigenvalue Study
%
% Compare:
%
%   h+/h- = 0.95
%   h+/h- = 0.70
%   h+/h- = 0.50
%
% =========================================================================

clear
clc
close all

%% Logging

problem_name = 'methodB_accuracy_eigenvalues';

fid = logger_init(problem_name);

logger(fid,'INFO',...
    'Starting Method B study');

%% Parameters

L = 1.0;
U = 1.0;
k = 2e-3;

bc.W.type  = 'dirichlet';
bc.W.value = 0.0;

bc.E.type  = 'dirichlet';
bc.E.value = 1.0;

ratio_list = [0.95 0.70 0.50];

%% Fixed mesh size

N = 50;

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

%% Storage

errors  = zeros(length(ratio_list),1);
minReal = zeros(length(ratio_list),1);

colors = lines(length(ratio_list));

%% =======================================================================
% Solution comparison figure
%% =======================================================================

fig_sol = figure('Color','white');
hold on

%% =======================================================================
% Eigenvalue comparison figure
%% =======================================================================

fig_eigs = figure('Color','white');
hold on

%% =======================================================================
% Loop over stretching ratios
%% =======================================================================

for ir = 1:length(ratio_list)

    r = ratio_list(ir);

    logger(fid,'INFO',...
        '================================================');

    logger(fid,'INFO',...
        'Stretch ratio h+/h- = %.2f',r);

    %% ---------------------------------------------------------
    % Mesh
    %% ---------------------------------------------------------

    dx0 = L*(1-r)/(1-r^N);

    dx = dx0*r.^(0:N-1);

    x = [0;cumsum(dx(:))];

    mesh.x = x;
    mesh.dx_w = [dx(1);dx(:)];
    mesh.dx_e = [dx(:);dx(end)];

    %% ---------------------------------------------------------
    % Numerical solution
    %% ---------------------------------------------------------

    [phi,A,rhs] = convection_diffusion_1d( ...
        mesh,...
        U,...
        k,...
        bc,...
        @convection_methodB_1d,...
        fid);

    %% ---------------------------------------------------------
    % Analytical solution (Eq. 1.3)
    %% ---------------------------------------------------------

    phi_exact = ...
        (exp(U*x/k)-1) ...
       /(exp(U*L/k)-1);

    %% ---------------------------------------------------------
    % Error
    %% ---------------------------------------------------------

    errors(ir) = ...
        sqrt(mean((phi-phi_exact).^2));

    logger(fid,'INFO',...
        'L2 error = %.6e',...
        errors(ir));

    %% ---------------------------------------------------------
    % Matrix assembly
    %% ---------------------------------------------------------

    n_internal = length(x)-2;

    M = zeros(n_internal,n_internal);

    for I = 1:n_internal

        i = I+1;

        [aW,aP,aE,~] = ...
            convection_methodB_1d( ...
            mesh,U,k,i);

        if I > 1
            M(I,I-1) = aW;
        end

        M(I,I) = aP;

        if I < n_internal
            M(I,I+1) = aE;
        end

    end

    eigvals = eig(M);

    minReal(ir) = ...
        min(real(eigvals));

    logger(fid,'INFO',...
        'Minimum Re(lambda) = %.6e',...
        minReal(ir));

    %% ---------------------------------------------------------
    % Solution comparison
    %% ---------------------------------------------------------

    figure(fig_sol)

    plot(x,...
         phi,...
         '-',...
         'LineWidth',2,...
         'Color',colors(ir,:));

    plot(x,...
         phi_exact,...
         '--',...
         'LineWidth',2,...
         'Color',colors(ir,:));

    %% ---------------------------------------------------------
    % Eigenvalues
    %% ---------------------------------------------------------

    figure(fig_eigs)

    plot(real(eigvals),...
         imag(eigvals),...
         'o',...
         'MarkerSize',6,...
         'LineWidth',1.5,...
         'Color',colors(ir,:));

    %% ---------------------------------------------------------
    % Save individual data
    %% ---------------------------------------------------------

    save( ...
        fullfile( ...
        data_dir,...
        sprintf('ratio_%03d.mat', ...
        round(100*r))), ...
        'r',...
        'phi',...
        'phi_exact',...
        'eigvals');

end

%% =======================================================================
% Finalize solution figure
%% =======================================================================

figure(fig_sol)

xlabel('x');
ylabel('\phi');

title('Method B: Numerical vs Analytical Solution');

legend( ...
    'Num. h+/h-=0.95',...
    'Exact h+/h-=0.95',...
    'Num. h+/h-=0.70',...
    'Exact h+/h-=0.70',...
    'Num. h+/h-=0.50',...
    'Exact h+/h-=0.50',...
    'Location','best');

grid on
box on

saveas( ...
    fig_sol,...
    fullfile( ...
    figures_dir,...
    'methodB_solution_comparison.png'));

%% =======================================================================
% Finalize eigenvalue figure
%% =======================================================================

figure(fig_eigs)

xline(0,'k--');

xlabel('Re(\lambda)');
ylabel('Im(\lambda)');

title('Method B Eigenvalue Spectra');

legend( ...
    'h+/h-=0.95',...
    'h+/h-=0.70',...
    'h+/h-=0.50',...
    'Location','best');

grid on
box on

saveas( ...
    fig_eigs,...
    fullfile( ...
    figures_dir,...
    'methodB_eigenvalue_comparison.png'));

%% =======================================================================
% Summary
%% =======================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' Method B Summary\n');
fprintf('====================================================\n');
fprintf('\n');
fprintf(' h+/h-      L2 Error      Min Re(lambda)\n');
fprintf('----------------------------------------\n');

for i=1:length(ratio_list)

    fprintf('%6.2f   %12.4e   %12.4e\n',...
        ratio_list(i),...
        errors(i),...
        minReal(i));

end

logger(fid,'INFO',...
    'Study completed successfully');

logger_close(fid,problem_name);
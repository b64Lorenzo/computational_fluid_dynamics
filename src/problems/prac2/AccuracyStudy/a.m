%% =========================================================================
% CFD - B3 Local Truncation Error Study
%
% Verifies that the B3 derivative approximation is second-order accurate.
%
% =========================================================================

clear
clc
close all

%% Logging

problem_name = 'b3_local_order';

fid = logger_init(problem_name);

logger(fid,'INFO','Starting B3 local order study');

%% Parameters

L = 1.0;
U = 1.0;
k = 1.0;

N_list = [20 40 80 160 320 640];

nCases = length(N_list);

h_vec = zeros(nCases,1);
err   = zeros(nCases,1);

%% Exact solution

for icase = 1:nCases

    N = N_list(icase);

    logger(fid,'INFO', ...
        'Running N=%d',N);

    h = L/N;

    h_vec(icase) = h;

    x = linspace(0,L,N+1)';

    % Exact solution Eq. (1.3)

    phi = ...
        (exp(U*x/k)-1) ...
        /(exp(U*L/k)-1);

    % Exact derivative

    dphi_exact = ...
        (U/k)*exp(U*x/k) ...
        /(exp(U*L/k)-1);

    % Use only points where B3 stencil exists

    idx = 3:N;

    dphi_b3 = ...
        (3*phi(idx) ...
        -4*phi(idx-1) ...
        +phi(idx-2)) ...
        /(2*h);

    err(icase) = sqrt(mean( ...
        (dphi_b3-dphi_exact(idx)).^2));

    logger(fid,'INFO', ...
        'h     = %.6e',h);

    logger(fid,'INFO', ...
        'Error = %.6e',err(icase));

end

%% Rates

rate = zeros(nCases-1,1);

for i=1:nCases-1

    rate(i) = ...
        log(err(i)/err(i+1)) ...
        / log(2);

end

%% Print table

fprintf('\n');
fprintf('===============================================\n');
fprintf(' B3 LOCAL ACCURACY STUDY\n');
fprintf('===============================================\n');
fprintf('\n');
fprintf('N        Error          Rate\n');
fprintf('-----------------------------------------------\n');

for i=1:nCases

    if i==1

        fprintf('%5d   %12.4e      ---\n', ...
            N_list(i),err(i));

    else

        fprintf('%5d   %12.4e   %8.4f\n', ...
            N_list(i),err(i),rate(i-1));

    end

end

%% Average order

p = mean(rate(end-2:end));

fprintf('\nObserved order = %.6f\n',p);

logger(fid,'INFO', ...
    'Observed order = %.6f',p);

%% Plot

figure('Color','white');

loglog( ...
    h_vec,...
    err,...
    'o-',...
    'LineWidth',2,...
    'MarkerSize',8);

hold on

ref = err(1) ...
      *(h_vec/h_vec(1)).^2;

loglog( ...
    h_vec,...
    ref,...
    'k--',...
    'LineWidth',2);

xlabel('h');
ylabel('Error');

title('B3 Local Truncation Error');

legend( ...
    'Numerical',...
    'O(h^2)',...
    'Location','best');

grid on
box on

%% Save

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

saveas( ...
    gcf,...
    fullfile(figures_dir,...
    'b3_local_accuracy.png'));

save( ...
    fullfile(data_dir,...
    'b3_local_accuracy.mat'), ...
    'N_list',...
    'h_vec',...
    'err',...
    'rate');

logger(fid,'INFO', ...
    'Study completed');

logger_close(fid,problem_name);
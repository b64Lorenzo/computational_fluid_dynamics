%% =========================================================================
% EXERCISE 3 - WAVE SIMULATION
%
% d(phi)/dt + U*d(phi)/dx = 0
%
% Explicit Euler + Upwind
%
% Exact solution:
%
% phi(x,t) = 2 + sin(6*pi*(x-U*t))
%
% (a) Verify stability condition:
%
%         eta <= 1
%
% (b) Investigate numerical damping
%
%% =========================================================================

clear
clc
close all

%% =========================================================================
% Logging
%% =========================================================================

problem_name = 'exercise3_wave_simulation';

fid = logger_init(problem_name);

logger(fid,'INFO','Starting Exercise 3');

%% =========================================================================
% Directories
%% =========================================================================

results_dir = fullfile('results',problem_name);

figures_dir = fullfile(results_dir,'figures');

data_dir = fullfile(results_dir,'data');

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

N = 201;

t_end = 1.0;

%% =========================================================================
% Mesh
%% =========================================================================

x = linspace(0,L,N)';

dx = x(2)-x(1);

%% =========================================================================
% Initial condition
%% =========================================================================

phi0 = 2 + sin(6*pi*x);

amp0 = 0.5*(max(phi0)-min(phi0));

%% =========================================================================
% CFL Numbers
%% =========================================================================

eta_list = [ ...
    0.25 ...
    0.50 ...
    0.75 ...
    0.90 ...
    0.95 ...
    1.00 ...
    1.05 ...
    1.10 ...
    1.20 ];

nCases = length(eta_list);

final_error = zeros(nCases,1);

amplitude_ratio = zeros(nCases,1);

all_phi_final = cell(nCases,1);

%% =========================================================================
% Simulations
%% =========================================================================

for icase = 1:nCases

    eta = eta_list(icase);

    dt = eta*dx/U;

    time = 0:dt:t_end;

    nt = length(time);

    phi = zeros(N,nt);

    phi(:,1) = phi0;

    for n = 1:nt-1

        tnp1 = time(n+1);

        %% Boundary conditions

        phi(1,n+1) = ...
            2 + sin(6*pi*(0-U*tnp1));

        phi(end,n+1) = ...
            2 + sin(6*pi*(1-U*tnp1));

        %% Explicit Upwind

        for i = 2:N-1

            phi(i,n+1) = ...
                phi(i,n) ...
              - eta * ...
                (phi(i,n)-phi(i-1,n));

        end

    end

    %% Exact solution

    phi_exact = ...
        2 + sin(6*pi*(x-U*t_end));

    %% Error

    final_error(icase) = ...
        sqrt(mean( ...
        (phi(:,end)-phi_exact).^2));

    %% Numerical damping

    amp_final = ...
        0.5*(max(phi(:,end))-min(phi(:,end)));

    amplitude_ratio(icase) = ...
        amp_final/amp0;

    all_phi_final{icase} = phi(:,end);

end

%% =========================================================================
% Figure 1
% Stable solutions only
%% =========================================================================

stable_idx = eta_list <= 1.0;

phi_exact = ...
    2 + sin(6*pi*(x-U*t_end));

figure('Color','w')

plot(x,...
     phi_exact,...
     'k--',...
     'LineWidth',3)

hold on

cols = lines(sum(stable_idx));

counter = 1;

for icase = 1:nCases

    if stable_idx(icase)

        plot(x,...
             all_phi_final{icase},...
             'Color',cols(counter,:),...
             'LineWidth',2)

        counter = counter + 1;

    end

end

xlabel('x')
ylabel('\phi')

title('Stable Solutions (\eta \leq 1)')

ylim([0.8 3.2])

grid on
box on

legend_entries = {'Exact'};

for icase = find(stable_idx)

    legend_entries{end+1} = ...
        sprintf('\\eta = %.2f',eta_list(icase));

end

legend(legend_entries,...
       'Location','best')

saveas(gcf,...
       fullfile(figures_dir,...
       'stable_wave_profiles.png'))

%% =========================================================================
% Figure 2
% Stability and damping
%% =========================================================================

figure('Color','w')

subplot(1,2,1)

semilogy(eta_list,...
         final_error,...
         'o-',...
         'LineWidth',2,...
         'MarkerSize',8)

hold on

xline(1,...
      'r--',...
      'LineWidth',2)

xlabel('\eta')
ylabel('L_2 Error')

title('Stability Investigation')

grid on
box on

subplot(1,2,2)

stable_amplitude = ...
    amplitude_ratio(eta_list<=1);

stable_eta = ...
    eta_list(eta_list<=1);

plot(stable_eta,...
     stable_amplitude,...
     'o-',...
     'LineWidth',2,...
     'MarkerSize',8)

hold on

yline(1,...
      'k--',...
      'LineWidth',2)

xlabel('\eta')
ylabel('Amplitude Ratio')

title('Numerical Damping')

grid on
box on

saveas(gcf,...
       fullfile(figures_dir,...
       'stability_and_damping.png'))

%% =========================================================================
% Results
%% =========================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('EXERCISE 3 RESULTS\n');
fprintf('====================================================\n');

for icase = 1:nCases

    fprintf( ...
        'eta = %.2f | Error = %.6e | AmpRatio = %.6f\n',...
        eta_list(icase),...
        final_error(icase),...
        amplitude_ratio(icase));

end

fprintf('====================================================\n');

fprintf('\n');
fprintf('Part (a):\n');
fprintf('Stable for eta <= 1.\n');
fprintf('Unstable for eta > 1.\n');

fprintf('\n');
fprintf('Part (b):\n');
fprintf('Minimum damping obtained near eta = 1.\n');

%% =========================================================================
% Save Data
%% =========================================================================

save(fullfile(data_dir,...
     'exercise3_results.mat'),...
     'eta_list',...
     'final_error',...
     'amplitude_ratio',...
     'all_phi_final');

logger(fid,'INFO','Exercise completed');

logger_close(fid,problem_name);
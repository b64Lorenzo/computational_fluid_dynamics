%% =========================================================================
% Computational Fluid Dynamics
%
% Project startup script
% =========================================================================

clearvars

project_root = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(project_root,'src')));

disp('===========================================================');
disp(' Computational Fluid Dynamics Project Initialized');
disp('===========================================================');
disp(['Project Root: ', project_root]);
disp('Source directories added to MATLAB path.');
disp('===========================================================');
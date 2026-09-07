clear
clc

problem_name = 'logger_test';

fid = logger_init(problem_name);

logger(fid,'INFO','L = %.2f',10);
logger(fid,'INFO','U = %.2f',50);
logger(fid,'INFO','k = %.3e',1e-3);

logger(fid,'INFO','Creating mesh');
logger(fid,'INFO','Assembling matrix');
logger(fid,'INFO','Solving system');
logger(fid,'INFO','Saving results');

logger_close(fid,problem_name);
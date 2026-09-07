function fid = logger_init(problem_name)

timestamp = datestr(now,'yyyymmdd_HHMMSS');

log_dir = fullfile('results',problem_name,'logs');

if ~exist(log_dir,'dir')
    mkdir(log_dir);
end

log_file = fullfile( ...
    log_dir, ...
    sprintf('%s_%s.log',problem_name,timestamp));

fid = fopen(log_file,'w');

logger(fid,'INFO', ...
    '============================================================');

logger(fid,'INFO', ...
    'STARTING %s',upper(problem_name));

logger(fid,'INFO', ...
    '============================================================');

logger(fid,'INFO','Log file: %s',log_file);

end
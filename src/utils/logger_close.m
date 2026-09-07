function logger_close(fid,problem_name)

logger(fid,'INFO',...
    '================================================================================');

logger(fid,'INFO',...
    '%s FINISHED', ...
    upper(problem_name));

logger(fid,'INFO',...
    '================================================================================');

fclose(fid);

end

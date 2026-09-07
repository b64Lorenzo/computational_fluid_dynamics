function logger_close(fid,problem_name)

logger(fid,'INFO', ...
    '============================================================');

logger(fid,'INFO', ...
    '%s FINISHED', ...
    upper(problem_name));

logger(fid,'INFO', ...
    '============================================================');

if ~isempty(fid) && fid > 0
    fclose(fid);
end

end
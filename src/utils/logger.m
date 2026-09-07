function logger(fid,level,msg,varargin)

message = sprintf(msg,varargin{:});

timestamp = datestr(now,'yyyy-mm-dd HH:MM:SS,FFF');

line = sprintf('%s | %-7s | %s', ...
    timestamp,...
    upper(level),...
    message);

fprintf('%s\n',line);

if ~isempty(fid) && fid > 0
    fprintf(fid,'%s\n',line);
end

end
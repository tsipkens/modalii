
% WRITE_JSON  Write JSON structured configuration files. 
%  
%  AUTHOR: Timothy Sipkens, 2023-10-12

function j = write_json(prop, file)

if ~exist('file', 'var'); file = []; end

% Remove options, as not useful universally. 
prop = rmfield(prop, 'opts');

f = fieldnames(prop);

% Convert function handles to strings for saving. 
for ii=1:length(f)
    if isa(prop.(f{ii}), 'function_handle')
        prop.(f{ii}) = func2str(prop.(f{ii}));
    end
end

% JSON encode data. 
j = jsonencode(prop);

% Write file IF a file name (i.e., FILE) was provided. 
if ~isempty(file)
    fid = fopen(file, 'w');
    fprintf(fid, '%s', j);
    fclose(fid);
end

end


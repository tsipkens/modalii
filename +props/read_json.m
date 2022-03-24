
% READ_JSON  Read JSON structured configuration files. 
%  Allows for C++ or Javscript style commenting.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-20

function s = read_json(file)

fid = fopen(file);
raw = fread(fid, inf);  % raw file contents
str = char(raw');  % transform to char
fclose(fid);

% Remove comments.
str = erase(erase(eraseBetween( ...
    erase(eraseBetween(str, "//", newline), "//"), ...
    "/*", "*/"), "/*"), "*/");

s = jsondecode(str);

% Attempt to interpret Matlab expressions.
fiel = fields(s);
for ii=1:length(fiel)
    t0 = s.(fiel{ii});
    
    if isa(t0, 'char')
        [converted, success] = str2num(t0);
        if success
            s.(fiel{ii}) = converted;
        end
    end
end

end


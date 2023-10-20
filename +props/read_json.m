
% READ_JSON  Read JSON structured configuration files. 
%  Allows for C++ or Javscript style commenting.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-20

function prop = read_json(file)

fid = fopen(file);
raw = fread(fid, inf);  % raw file contents
str = char(raw');  % transform to char
fclose(fid);

% Remove comments.
str = erase(erase(eraseBetween( ...
    erase(eraseBetween(str, "//", newline), "//"), ...
    "/*", "*/"), "/*"), "*/");

prop = jsondecode(str);  % use "prop" as that is what is used in JSON

% Attempt to interpret Matlab expressions.
fiel = fields(prop);
for ii=1:length(fiel)
    t0 = prop.(fiel{ii});
    
    if isa(t0, 'char')
        if strcmp(t0(1), '@')  % if function handle then evaluate
            prop.(fiel{ii}) = eval(t0);

        else  % otherwise, try to convert to number
            [converted, success] = str2num(t0);
            if success
                prop.(fiel{ii}) = converted;
            end
        end
    end
end

end


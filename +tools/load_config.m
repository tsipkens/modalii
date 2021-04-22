
% LOAD_CONFIG  Loads settings from configuration file (YML, YAML, or JSON). 
%  Files are loaded in order supplied, overwriting properties where
%  relevant. 
%  
%  AUTHOR: Timothy Sipkens, 2021-03-25

function prop = load_config(fnames, prop)

if ~iscell(fnames); fnames = {fnames}; end

if ~exist('prop', 'var'); prop = []; end
if isempty(prop); prop = struct(); end

for ii=1:length(fnames)
    
    prop0 = read_json(fnames{ii});  % read new settings
    
    
    f = fieldnames(prop0);
    
    % Copy (or overwrite) existing settings.
    for jj = 1:length(f)
        
        % Attempt to interpret Matlab expressions.
        prop0.(f{jj}) = interpret(prop0.(f{jj}));
        
        % Overwrite existing configuration.
        prop.(f{jj}) = prop0.(f{jj});
    end
    
    %-- Handle if/then/else statements ------%
    idx = find(contains(f, "if"));
    
    % Look for if... then... else pattern.
    idx(idx == length(f)) = [];  % remove from flag if last item
    idx(~contains(f(idx + 1), "then")) = [];
    idx(~contains(f(idx + 2), "else")) = [];
    
    for ii=idx
        if prop0.(f{ii})
            fj = fieldnames(prop0.(f{ii+1}));
            for jj=1:length(fj)
                t0 = interpret(prop0.(f{ii+1}).(fj{jj}));
                prop.(fj{jj}) = t0;
            end
        else
            fj = fieldnames(prop0.(f{ii+2}));
            for jj=1:length(fj)
                t0 = interpret(prop0.(f{ii+2}).(fj{jj}));
                prop.(fj{jj}) = t0;
            end
        end
        
        prop = rmfield(prop, f{ii});
        prop = rmfield(prop, f{ii+1});
        prop = rmfield(prop, f{ii+2});
    end
    
end

end



% INTERPRET  Attempt to interpret Matlab expressions.
function e0 = interpret(f)

if isa(f, 'char')
    [e0, success] = str2num(f);
    if success; return; end

    try  % try 'eval(...)'
        e0 = eval(f);
    catch
        e0 = f;
        return;  % just continue
    end
else
    e0 = f;
end

end




% READ_JSON  Read JSON structured configuration files. 
%  Allows for C++ or Javscript style commenting.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-20

function results = read_json(file)

fid = fopen(file);
raw = fread(fid, inf);  % raw file contents
str = char(raw');  % transform to char
fclose(fid);

% Remove comments.
str = erase(erase(eraseBetween( ...
    erase(eraseBetween(str, "//", newline), "//"), ...
    "/*", "*/"), "/*"), "*/");

results = jsondecode(str);

end




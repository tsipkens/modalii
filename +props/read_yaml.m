
% READ_YAML  Prase YAML files into Matlab structures. 
%  
%  See also YAML reader by Serge, which informed on this code. 
%  https://github.com/serg3y/MatLab-YAML.
%  
%  AUTHOR: Timothy Sipkens, 2023-10-18

function [prop] = read_yaml(txt, prop)

if ~exist('prop', 'var'); prop = []; end
if isempty(prop); prop = struct(); end

% Convert YAML text to data (uses SnakeYAML).
if ~any(contains(javaclasspath('-all'),'snakeyaml'))
    yamlsetup;
end

% Read YAML as text, if not already YAML formatted text. 
% Allows for input as YAML text. 
if isfile(txt)
    txt = fileread(txt);
end

% Parse YAML text.
J = org.yaml.snakeyaml.Yaml().load(txt);  % load as java class
prop_new = java2matlab(J, 1, 0);  % convert to MatLab (recursive, general parse)

% Update existing prop structure, if supplied.
% This is necessary to use pre-existing quantities in parsing function handles. 
fields = fieldnames(prop_new);
for ii=1:length(fields)
    prop.(fields{ii}) = prop_new.(fields{ii});
end

prop = customparse(prop); % custom parsing specific to this format (recursive)

end


% CUSTOMPARSE =========================================%
% A specific YAML parser to further interpret data (e.g., evaluating
% function handles stored in results). 
function prop = customparse(prop)

% Attempt to further interpret Matlab expressions.
fiel = fields(prop);
for ii=1:length(fiel)
    t0 = prop.(fiel{ii});

    if isa(t0, 'struct')
        for jj=1:length(t0)  % loop through struct (e.g., for evap model)
            prop.(fiel{ii})(jj) = customparse(t0(jj));
        end

    elseif isa(t0, 'char')
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


% JAVA2MATLAB =========================================%
% When JoinCells=true nested cells are joined to form an array with
% any number of dimension. During recursive call through the nested cells
% 'depth' is used to track recursion depth and hence dimension number,
% 'maxdepth' tracks total number of dimensions. Before exiting from
% depth==1 permute is used to swap dimension order.
function [data,maxdepth] = java2matlab(J, join, depth)

maxdepth = depth; %init
switch class(J)
    case 'java.util.ArrayList' %java list
        if J.size == 0
            data = {};  % convert to empty cell

        else
            data = J.toArray.cell';  % convert to vector of cells
            if join
                types = {'double' 'logical' 'java.util.ArrayList' 'java.util.LinkedHashMap'}; %data types that can be joined
                ValidType = ismember(class(data{1}),types); %is first element of a type that can be joined?
                SameType = all(cellfun(@(x)isequal(class(x),class(data{1})),data)); %are all elements of the same type?
                SameSize = all(cellfun(@(x)isequal(size(x),size(data{1})),data)); %are all elements of the same length?
                depth = (depth+1) * (ValidType && SameType && SameSize); %reset depth if joining is not possible, else elements will be joined
            end

            for k = 1:numel(data)
                [data{k},maxdepth] = java2matlab(data{k}, join, depth);  % receptively convert each element
            end

            if join && ValidType && SameType && SameSize
                try  % cell contents may have changed after recursive joining, field names may not match
                    data = cat(depth, data{:});  % join
                end
                if depth==1
                    if maxdepth==1
                        ord = [2 1];
                    else
                        ord = max(maxdepth,2):-1:1;
                        ord = ord([2 1 3:end]);
                    end
                    data = permute(data, ord);
                end
            end
        end

    case 'java.util.LinkedHashMap'
        val = J.values.toArray; %values
        if isempty(val)
            data = struct();
        else
            par = J.keySet.toArray.string; %slow but handles edge cases such as 'a: 1\n2: 2'
            % par = cell(1,J.size); t = J.keySet.iterator; for k = 1:J.size, par{k} = t.next.string; end %fast but fails on edge cases
            %the alternative is slower: par = cellstr(char(J.keySet.toArray)); %medium but fails on some edge cases
            par = matlab.lang.makeValidName(par); %ensure field names are valid, replace bad characters with _, prefix x to leading numbers
            par = matlab.lang.makeUniqueStrings(par); %ensure parameters are unique, append _1 _2 if needed
            for k = 1:numel(par)
                data.(par{k}) = java2matlab(val(k), join, depth); %assign and also convert contents
            end
        end

    case 'java.util.Date'
        data = datetime(J.getTime/1000, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
        data.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''';  % display format only, note java rounds milliseconds
        
    case {'char' 'double' 'logical'}
        data = J;  % do nothing, data is directly available
        
    otherwise
        data = J;
        fprintf(2,'Unsupported data type: %s\n',class(J));
end
end


% YAMLSETUP =========================================%
% Download and add SnakeYAML.jar to java class path.
% Download jar to same folder as this m file.
function yamlsetup()

fold = fileparts(mfilename('fullpath'));

jar = 'snakeyaml-2.0.jar';
pth = fullfile(fold,jar);
url = 'https://repo1.maven.org/maven2/org/yaml/snakeyaml/2.0/snakeyaml-2.0.jar';

% Download snakeyaml.
if ~isfile(pth)
    websave(pth,url);
end

% Add jar file temporarily to dynamic javaclasspaths.
if ~any(contains(javaclasspath('-all'), jar))
    javaaddpath(pth);
end
end


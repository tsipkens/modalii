
% WRITE_YAML  Prase YAML files into Matlab structures. 
%  
%  See also YAML reader by Serge, which informed on this code. 
%  https://github.com/serg3y/MatLab-YAML.
%  
%  See also TAML dump by Martin Koch.
%  https://github.com/MartinKoch123/yaml
%  
%  AUTHOR: Timothy Sipkens, 2023-10-18

function [txt] = write_yaml(prop, file)

if ~exist('file', 'var'); file = []; end

% Convert YAML text to data (uses SnakeYAML).
if ~any(contains(javaclasspath('-all'),'snakeyaml'))
    yamlsetup;
end

javaData = convert(prop);

dumperOptions = org.yaml.snakeyaml.DumperOptions();
txt = char(org.yaml.snakeyaml.Yaml(dumperOptions).dump(javaData));
txt(end) = []; % remove trailing line break

txt = replace(txt, '$%&?', 'null');
txt = replace(txt, "'", '"');

% Write file IF a file name (i.e., FILE) was provided. 
if ~isempty(file)
    fid = fopen(file, 'w');
    fprintf(fid, '%s', txt);
    fclose(fid);
end

end


% SUPPLEMENTAL FUNCTIONS ==============================%
function result = convert(data)
    if iscell(data)
        result = convertCell(data);
    elseif ischar(data) && isvector(data)
        result = convertString(data);
    elseif ~isscalar(data)
        result = convertArray(data);
    elseif isstruct(data)
        result = convertStruct(data);
    elseif isfloat(data)
        result = java.lang.Double(data);
    elseif isa(data, "int64")
        result = java.lang.Long(data);
    elseif isa(data, "uint32") || isa(data, "uint64")
        hexStr = dec2hex(data);
        result = java.math.BigInteger(hexStr, 16);
    elseif isinteger(data)
        result = java.lang.Integer(data);
    elseif islogical(data)
        result = java.lang.Boolean(data);
    elseif isstring(data)
        result = convertString(data);
    elseif isa(data,'function_handle')
        result = convertString(func2str(data));
    else
        error("yaml:dump:TypeNotSupported", "Data type '%s' is not supported.", class(data))
    end
end

function result = convertString(data)
    if contains(data, "$%&?")
        error("yaml:dump:NullPlaceholderNotAllowed", "Strings must not contain '%s' since it is used as a placeholder for null values.", "$%&?")
    end
    result = java.lang.String(data);
end

function result = convertStruct(data)
    result = java.util.LinkedHashMap();
    for key = string(fieldnames(data))'
        value = convert(data.(key));
        result.put(key, value);
    end
end

function result = convertCell(data)
    data = nest(data);
    result = java.util.ArrayList();
    for i = 1:length(data)
        result.add(convert(data{i}));
    end
end

function result = convertArray(data)
    result = convertCell(num2cell(data));
end

function result = nest(data)
    if isvector(data) || isempty(data)
        result = data;
        return
    end
    n = size(data, 1);
    nDimensions = length(size(data));
    result = cell(1, n);
    if nDimensions == 2
        for i = 1:n
            result{i} = data(i, :);
        end
    elseif nDimensions == 3
        for i = 1:n
            result{i} = squeeze(data(i, :, :));
        end
    else
        error("yaml:dump:HigherDimensionsNotSupported", "Arrays with more than three dimensions are not supported. Use nested cells instead.")
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


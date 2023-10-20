
% GET  Function to search for and obtain material properties.
%   
%  AUTHOR: Timothy Sipkens, 2023-10-17

function prop = get(spec, prop, opts)

if ~exist('opts', 'var'); opts = []; end  % used in MATLAB file load calls
if ~exist('prop', 'var'); prop = []; end
if isempty(prop); prop = props.init(); end

if ~iscell(spec); spec = {spec}; end

for jj=1:length(spec)
    if isfile(['+props/yaml/', spec{jj}, '.yaml'])  % YAML first
        prop = props.read_yaml(['+props/yaml/', spec{jj}, '.yaml'], prop);
        % NOTE: Property overriding occurs in read_yaml().
    
    elseif isfile(['+props/json/', spec{jj}, '.json'])  % JSON second
        prop0 = props.read_json(['+props/json/', spec{jj}, '.json']);
        
        % Update existing prop structure. 
        fields = fieldnames(prop0);
        for ii=1:length(fields)
            prop.(fields{ii}) = prop0.(fields{ii});
        end
    
    elseif isfile(['+props/', spec{jj}, '.m'])  % then MATLAB file
        prop = eval(['props.', spec{jj}, '(prop, opts);']);
    
    else
        error('Properties not available.');
    
    end
end

% If evaporation/sublimation model needs to be built.
% Relevant to multi-species evaporation.
if isfield(prop.opts, 'build_evap')
    if prop.opts.build_evap
        prop = props.build_evap(prop);
    end
end

end

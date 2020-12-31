
% PARSE_VARARGIN  Parses varargin inputs to class definitions and other functions.
% AUTHOR: Timothy Sipkens, 2019-06-25
%=========================================================================%

function [instance] = parse_varargin(instance, varargin)

ii = 1;
while ii<=length(varargin) % incorporate list of properties
    if isprop(instance, varargin{ii}) % manually set property
        instance.(varargin{ii}) = varargin{ii+1};
        ii = ii+2; % skip an input

    else  % incorporate opts given as struct
        aa = fieldnames(varargin{ii});
        bb = varargin{ii};
        for jj = 1:length(aa)
            if isfield(instance.opts,aa{jj})
                instance.opts.(aa{jj}) = bb.(aa{jj});
            end
        end
        ii = ii+1;
    end
end

end


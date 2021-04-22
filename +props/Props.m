
% PROPS  A container for for material properties. 
%  Allows for changes to function_handle parameters, which is required
%  for model selection of functions.
%  
%  PROP = Props(P) transforms input stuct P into an equivalent Props
%  object, which allows for changing of function_handle sub-parameters.
%  
%  ------------------------------------------------------------------------
%  
%  NOTE: Use of this class is substantially slower than using simple
%  builtin Matlab structs. 
%  
%  NOTE: Connections between other variables in class are established using
%  "prop.*" notation. For example, "prop.phi" is used to refer to the "phi"
%  field.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-21


classdef Props < dynamicprops
    
    properties
        % None, to start. 
        % Added in constructor method.
    end
    
    methods
        %== PROPS ========================================================%
        function prop = Props(p)
        f = fieldnames(p);

        for ii=1:length(f)
            addprop(prop, f{ii});

            if isa(p.(f{ii}), 'function_handle')
                prop.(f{ii}) = eval(func2str(p.(f{ii})));
            else
                prop.(f{ii}) = p.(f{ii});
            end
        end
        end
        
        %== CONVERT2STRUCT ===============================================%
        function p = convert2struct(prop)
        % Converts a Props object back to a simple Matlab struct. 
        
        f = properties(prop);
        for ii=1:length(f)
            p.(f{ii}) = prop.(f{ii});
        end
            
        end
    end
end


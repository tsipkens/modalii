
% BUILD_LIKELIHOOD  Build likelihood based on model and data inputs.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-20

function like = build_like(A, b, Lb, minimizer)

%-- Parse standard deviations ---------%
if ~exist('Lb', 'var'); Lb = []; end
if isempty(Lb); Lb = speye(length(b)); end

if ~exist('minimizer', 'var'); minimizer = []; end
if isempty(minimizer); minimizer = 'default'; end


if ~isa(b, 'function_handle')
    switch minimizer
        case {'default','levenberg-marquardt','vector','vector-xpr'} % minimization considers vector
            like = @(x) ((Lb')*(A(x)-b))./2;

        case {'fminsearch','scalar'} % minimization considers scalar
            like = @(x) norm(((Lb')*(A(x)-b))./2,2);

        otherwise
            disp('ERROR: Invalid minimization technique!');
            return;
    end
    
else
    if ~isa(Lb, 'function_handle')
        Lb = @(x) Lb;
    end
    switch minimizer
        case {'default','levenberg-marquardt','vector','vector-xpr'} % minimization considers vector
            like = @(x) ((Lb(x)')*(A(x)-b(x)))./2;

        case {'fminsearch','scalar'} % minimization considers scalar
            like = @(x)norm(((Lb(x)')*(A(x)-b(x)))./2,2);
            
        otherwise
            disp('ERROR: Invalid minimization technique!');
            return;
    end
end

end


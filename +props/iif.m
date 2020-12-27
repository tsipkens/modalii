
% IIF  If function for writing inline conditional statements.
% Author: Timothy Sipkens, 2020-12-27
%=========================================================================%

function out = iif(cond, a, b)
    out = b;
    out(cond) = a(cond);
end

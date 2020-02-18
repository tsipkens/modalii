
function [] = N2O(prop,~)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('prop','var'); prop = struct(); end
%-------------------------------------------------------------------------%

prop.mg = 7.30637e-26;
prop.gamma1 = 7/5;
prop.gamma2 = @(T)(prop.gamma1+1)/(prop.gamma1-1);

end



function prop = Ne(prop,~)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('prop','var'); prop = struct(); end
%-------------------------------------------------------------------------%

prop.mg = 3.35082e-26;
prop.gamma1 = 5/3;
prop.gamma2 = @(T)(prop.gamma1+1)/(prop.gamma1-1);

end


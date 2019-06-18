function [Gx,R,s] = credLinear(stats,jcb)

Gx = inv(((jcb')*jcb));
s = sqrt(diag(Gx));
R = bsxfun(@rdivide,Gx,s);
R = bsxfun(@rdivide,R,s');
R = full(R); % posterior correlation
s = full(s); % posterior standard deviation

end


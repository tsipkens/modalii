function [out] = plotProp(prop,propName,opts)

if nargin<3;
    opts.limits = [1000,3000];
elseif not(isfield(opts,'limits'))
    opts.limits = [1000,3000];
end

[outX,outY] = fplot(prop.(propName),opts.limits);
if not(isfield(opts,'plot'));
    plot(outX,outY);
elseif strcmp(opts.plot,'semilogy');
    semilogy(outX,outY);
elseif strcmp(opts.plot,'semilogx');
    semilogx(outX,outY);
else
    disp('Error plotting...');
end
out = [outX,outY];

end


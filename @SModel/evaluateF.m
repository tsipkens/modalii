function [Jout,mp] = evaluateF(smodel,x)
% EVALAUTEF Evaluate spectroscopic forward model (T/htmodel -> J). 

htmodel = smodel.htmodel; % embedded heat transfer model
prop = smodel.prop; % material properties

if nargin > 1 % update x values
    if length(x)==length(smodel.x)
        for ii=1:length(smodel.x)
            prop.(smodel.x{ii}) = x(ii);
        end
    else
        disp('Warning: QoIs parameter size mismatch.');
    end
end


%-- Consider size distribution -------------------------------------------%
if prop.sigma > 0.005 % currently models with lognormal
    dp0 = prop.dp0; % store current dp0 and apply as geometric mean
    
    %-- Discretize the space -----------------------%
    %   To do: Update to proper "logspace" discretization?
    dp_1 = logninv(0.01,real(log(prop.dp0)),prop.sigma);
    dp_2 = logninv(0.9999,real(log(prop.dp0)),prop.sigma);
    dp_x_diff = (dp_2-dp_1)/50;
    dp_x = (dp_1:dp_x_diff:dp_2)'; % differentiated space
    
    %-- Evaluate temperature and incandescence -------%
    T = htmodel.de_solve(dp_x);
    Jout = bsxfun(@times,(lognpdf(dp_x,log(dp0),prop.sigma).*dp_x_diff.*...
        pi.*(dp_x.*1e-9).^3)',...
        smodel.FModel(T,prop.Em)); % evaluate forward model for J
    Jout = sum(Jout,2);
    Jout = Jout.*prop.C_J; % scale incandescence for stability
    mp = []; % not currently calculating mass ratio for distribution
    
    
%-- Simple evaluation for monodisperse case ------------------------------%
else % for monodisperse case, simply evaluate the ODE directly
    [T,~,mp] = htmodel.de_solve(prop.dp0); % solve heat transfer model at dp0
    
    Jout = smodel.FModel(T,prop.Em); % evaluate forward model for J
    Jout = Jout.*prop.C_J; % scale incandescence by corresponding factor
    
    if strcmp(smodel.opts.multicolor,'constC-mass') % scale incandescence according to mass loss
        mpr = real(mp./mp(1));
        Jout = bsxfun(@times,Jout,mpr);
    end
    
end
%-------------------------------------------------------------------------%

end

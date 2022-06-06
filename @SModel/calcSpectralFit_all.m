function [To, Co, s_T,out] = calcSpectralFit_all(smodel,J)
% Spectral fitting, simultaneously

prop = smodel.prop;
l = smodel.l;
ntime = length(J(:,1,1)); % number of times
nlambda = length(J(1,1,:)); % number of wavelengths

s = std(J,[],2)./sqrt(length(J(1,:,1)));
J = mean(J,2);
data = reshape(J,[ntime*nlambda,1]);
nn = length(data);
data_std = reshape(s,[ntime*nlambda,1]);
data_Li = spdiags(ones(nn,1)./data_std,0,nn,nn);

opts.bFun = 0;
switch smodel.opts.multicolor
    case 'default'
        x0 = [prop.C_J.*ones(1,ntime),3000.*ones(1,ntime)];
        bb1 = @(x) reshape((bsxfun(@times,...
            bsxfun(@times,x(1:ntime)',smodel.blackbody(x(ntime+1:end)',l)),...
            prop.Em(l,prop.dp0)./(l.*1e-9)./smodel.data_sc)),[ntime*nlambda,1]);
            % scaled for stability, extra l from size parameter
        stats = Stats(bb1,[],opts);
    case 'constC'
        x0 = [prop.C_J,3000.*ones(1,ntime)];
        bb1 = @(x) reshape((bsxfun(@times,x(1).*smodel.blackbody(x(2:end)',l),...
            prop.Em(l,prop.dp0)./(l.*1e-9)./smodel.data_sc)),[ntime*nlambda,1]);
        stats = Stats(bb1,[],opts);
    case 'constC-mass'
        smodel.htmodel.Ti = 3000;
        [~,Mt] = smodel.htmodel.deSolve;% models loss in fv > reduced C_J
        Mt = real(Mt./Mt(1));
        x0 = [prop.C_J,3000.*(1-(1:ntime)./1000)];
        bb1 = @(x) reshape((bsxfun(@times,...
            bsxfun(@times,x(1).*Mt,smodel.blackbody(x(2:end)',l)),...
            prop.Em(l,prop.dp0)./(l.*1e-9)./smodel.data_sc)),[ntime*nlambda,1]);
        stats = Stats(bb1,[],opts);
    case 'priorC-smooth'
        x0 = [prop.C_J.*ones(1,ntime),3000.*ones(1,ntime)];
        bb1 = @(x) reshape((bsxfun(@times,...
            bsxfun(@times,x(1:ntime)',smodel.blackbody(x(ntime+1:end)',l)),...
            prop.Em(l,prop.dp0)./(l.*1e-9)./smodel.data_sc)),[ntime*nlambda,1]);
            % scaled for stability, extra l from size parameter
        x_pr = zeros(1,2*ntime);
        % x_pr = x0;
        alpha = 30;
        %Lx_i = spdiags(1./(0.3.*x0((ntime+1):end)),[zeros(1,ntime)]',0,2*ntime,2*ntime); % L2 regularization
        Lx_i = spdiags([ones(1,ntime),zeros(1,ntime)./2]',0,2*ntime,2*ntime); % 1st order Tikhonov
        Lx_i = Lx_i-spdiags([ones(1,ntime),zeros(1,ntime)./2]',1,2*ntime,2*ntime);
        Lx_i = Lx_i.*alpha;
        opts.minimize = 'vector-xpr';
        stats = Stats(bb1,[],opts,'Lx_ipr',Lx_i,'x_pr',x_pr);
end

stats.b = @(varargin)data;
stats.Lb_i = @(varargin)data_Li;
[mle,jcb] = stats.minimize(x0);
G_T = inv((jcb')*jcb);
out.resid = sqrt(sum(reshape(...
    (((stats.Lb_i(mle)')*(stats.model(mle)-stats.b(mle)))./2)...
    .^2,[ntime,nlambda]),2));

switch smodel.opts.multicolor
    case {'default','priorC-smooth'}
        To = mle((ntime+1):end);
        Co = mle(1:ntime);
        out.s_T = sqrt(diag(G_T((ntime+1):end,(ntime+1):end)));
        out.s_C = sqrt(diag(G_T(1:ntime,1:ntime)));
        R_TC = corrcov(G_T);
        out.r_TC = diag(R_TC,ntime+1);
    case {'constC','constC-mass'}
        To = mle(2:end);
        if strcmp('constC-mass',smodel.opts.multicolor); 
            Co = mle(1).*ones(ntime,1).*Mt;
        else
            Co = mle(1).*ones(ntime,1);
        end
        s_T = sqrt(diag(G_T(2:end,2:end)));
        out.s_C = sqrt(G_T(1,1));
        R_TC = corrcov(G_T);
        out.r_TC = R_TC(2:end,1);
end


end


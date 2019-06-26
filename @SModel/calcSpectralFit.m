function [Tout,Cout,s_T,out] = calcSpectralFit(smodel,J)
% Spectral fitting, sequential

ntime = length(J(:,1,1)); % number of times
s = std(J,[],2)./sqrt(length(J(1,:,1)));
J = mean(J,2);
nshots = length(J(1,:,1)); % number of shots

prop = smodel.prop;
l = smodel.l;
s_T = zeros(ntime,nshots);
s_C = zeros(ntime,nshots);

opts.bFun = 0;
switch smodel.opts.multicolor
    case 'default'
        x0 = [3000,prop.C_J];
        bb1 = @(x) x(2).*squeeze(smodel.blackbody(x(1),l))...
            .*(prop.Em(l,prop.dp0)./(l.*1e-9)./smodel.data_sc)';
            % scaled for stability, extra l from size parameter
        r_TC = zeros(ntime,nshots);
        stats = Stats(bb1,[],opts);
    case 'constC'
        x0 = 3000;
        bb1 = @(x) prop.C_J.*squeeze(smodel.blackbody(x(1),l))...
            .*(prop.Em(l,prop.dp0)./(l.*1e-9)./smodel.data_sc)';
        stats = Stats(bb1,[],opts);
    case 'constC-mass'
        x0 = 3000;
        smodel.htmodel.Ti = 2500;
        [~,Mt] = smodel.htmodel.deSolve;% models loss in fv > reduced C_J
        Mt = real(Mt./Mt(1));
        bb0 = @(x) prop.C_J.*squeeze(smodel.blackbody(x(1),l))...
            .*(prop.Em(l,prop.dp0)./(l.*1e-9)./smodel.data_sc)';
    case 'priorC'
        x0 = [3000,prop.C_J];
        bb1 = @(x) x(2).*squeeze(smodel.blackbody(x(1),l))...
            .*(prop.Em(l,prop.dp0)./(l.*1e-9)./smodel.data_sc)';
        r_TC = zeros(length(J(:,1,1)),length(J(1,:,1)));
        x_pr = [3000,prop.C_J];
        Lx_i = spdiags([0;1/(0.5.*prop.C_J)],0,...
            length(x0),length(x0));
        opts.minimize = 'vector-xpr';
        opts.prior = 'xpr';
        stats = Stats(bb1,[],opts,'Lx_ipr',Lx_i,'x_pr',x_pr);
    case 'vecC'
        x0 = 3000;
        bb0 = @(x) squeeze(smodel.blackbody(x(1),l))...
            .*(prop.Em(l,prop.dp0)./(l.*1e-9)./smodel.data_sc)';
        Mt = prop.C_J;
	case 'priorT'
        x0 = [3000,prop.C_J];
        bb1 = @(x) x(2).*squeeze(smodel.blackbody(x(1),l))...
            .*(prop.Em(l,prop.dp0)./(l.*1e-9)./smodel.data_sc)';
        r_TC = zeros(length(J(:,1,1)),length(J(1,:,1)));
        x_pr = [prop.Tg,prop.C_J];
        Lx_i = spdiags([1/100;0],0,...
            length(x0),length(x0));
        opts.minimize = 'vector-xpr';
        stats = Stats(bb1,[],opts,'Lx_ipr',Lx_i,'x_pr',x_pr);
    case 'priorC-nuisance'
        x0 = [3000,prop.C_J,prop.Em(l,prop.dp0)];
        bb1 = @(x) x(2).*squeeze(smodel.blackbody(x(1),l))...
            .*(x(3:end)./(l.*1e-9)./smodel.data_sc)';
        r_TC = zeros(length(J(:,1,1)),length(J(1,:,1)));
        x_pr = [3000,prop.C_J,prop.Em(l,prop.dp0)];
        Lx_i = spdiags([0;1/(0.3.*prop.C_J);1./(0.1.*prop.Em(l,prop.dp0))'],0,...
            length(x0),length(x0));
        opts.minimize = 'vector-xpr';
        opts.prior = 'xpr';
        stats = Stats(bb1,[],opts,'Lx_ipr',Lx_i,'x_pr',x_pr);
end

beta = zeros(ntime,nshots,length(x0));
resid = zeros(ntime,nshots);
for ii = 1:ntime % Time loop
    if or(strcmp(smodel.opts.multicolor,'constC-mass'),...
            strcmp(smodel.opts.multicolor,'vecC'))
        bb1 = @(x) bb0(x).*Mt(ii);
        stats = Stats(bb1,[],opts);
    end
    
    for jj = 1:nshots % Shot loop, only one if averaged
        data = squeeze(J(ii,jj,:));
        nn = length(data);
        data_std = squeeze(s(ii,1,:));
        data_std(data_std==0) = 1e-3; % remove zeros, may induce errors
        data_Li = spdiags(ones(nn,1)./data_std,0,nn,nn);
        
        stats.opts.bFun = 1;
        stats.b = @(varargin)data;
        stats.Lb_i = @(varargin)data_Li;
        [mle,jcb] = stats.minimize(x0);
        G_T = inv((jcb')*jcb);
        
        beta(ii,jj,:) = mle;
        s_T(ii,jj) = sqrt(G_T(1,1));
        
        if length(mle)==2
            r_TC(ii,jj) = G_T(1,2)./sqrt(G_T(1,1)*G_T(2,2));
            s_C(ii,jj) = sqrt(G_T(2,2));
        else
            s_C(ii,jj) = 0;
        end
        resid(ii,jj) = norm(((stats.Lb_i(mle)')*...
            (stats.model(mle)-stats.b(mle)))./2);
    end
    if ii==1
        fprintf('Calculating temperatures... ');
        tools.textbar(0);
    elseif (mod(ii,50)==0)
        tools.textbar(ii/length(J(:,1,1)));
        %disp([num2str(round(100*ii/length(J(:,1,1)),1)),'% complete.']);
    end
end

switch smodel.opts.multicolor
    case {'priorC','default','priorT'}
        Cout = beta(:,:,2);
        out.r_TC = r_TC;
        out.s_C = s_C;
        out.resid = resid;
    case {'constC'}
        out.r_TC = 0;
        out.s_C = s_C;
        out.resid = resid;
        Cout = prop.C_J.*ones(ntime,1);
    case {'constC-mass','vecC'}
        out.r_TC = 0;
        out.s_C = s_C;
        out.resid = resid;
        Cout = prop.C_J;
    case {'priorC-nuisance'}
        Cout = beta(:,:,2);
        out.r_TC = r_TC;
        out.s_C = s_C;
        out.resid = resid;
        out.Em = beta(:,:,3:end);
end
Tout = beta(:,:,1);

end



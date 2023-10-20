
% CALCSPECTRALFIT  Spectral fitting, sequential.
%  
%  AUTHOR: Timothy Sipkens, Refactoring started on 2023-10-20

function [To, Co, s_T, out] = calcSpectralFit(smodel, J)

ntime = length(J(:,1,1)); % number of times
s = std(J,[],2)./sqrt(length(J(1,:,1)));
J = mean(J,2);
nshots = length(J(1,:,1)); % number of shots

prop = smodel.prop;
if ~isfield(prop, 'C_J'); prop.C_J = []; end
if isempty(prop.C_J); prop.C_J = 1; end

l = smodel.l;
s_T = zeros(ntime,nshots);
s_C = zeros(ntime,nshots);

opts.bFun = 0;
switch smodel.opts.multicolor
    case 'default'
        x0 = [1e3, 0];
        bb1 = @(x) (10.^x(2)) .* ...
            squeeze(smodel.blackbody(x(1),l))...
            .*(prop.Em(l, prop.dp0, 1, prop)./(l.*1e-9))';
            % scaled for stability, extra l from size parameter
        r_TC = zeros(ntime,nshots);

    case 'constC'  % holds the pre-factor constant over the signal
        x0 = 1e3;
        bb1 = @(x) prop.C_J .* squeeze(smodel.blackbody(x(1),l))...
            .*(prop.Em(l, prop.dp0, 1, prop)./(l.*1e-9))';
end

beta = zeros(ntime,nshots,length(x0));
resid = zeros(ntime,nshots);
for ii = 1:ntime  % time loop
    
    for jj = 1:nshots % Shot loop, only one if averaged
        data = squeeze(J(ii,jj,:));
        nn = length(data);
        data_std = squeeze(s(ii,1,:));
        data_std(data_std==0) = 1e-3; % remove zeros, may induce errors
        data_Li = spdiags(ones(nn,1)./data_std,0,nn,nn);
        
        tmp = bb1(x0);
        like = @(x) (log(data ./ data(1)) - log(bb1(x) ./ tmp(1)))';
        
        [mle,jcb] = stats.minimize(x0, like);
        mle(2) = 10 .^ (mle(2) ./ tmp(1) .* data(1)); % restore constant to full values
        G_T = inv((jcb')*jcb);
        
        beta(ii,jj,:) = mle;
        s_T(ii,jj) = sqrt(G_T(1,1));
    end
    if ii==1
        fprintf('Calculating temperatures:');
        tools.textbar(0);
    elseif (mod(ii,50)==0)
        tools.textbar(ii/length(J(:,1,1)));
        %disp([num2str(round(100*ii/length(J(:,1,1)),1)),'% complete.']);
    end
end

switch smodel.opts.multicolor
    case {'priorC','default','priorT'}
        Co = beta(:,:,2);
    case {'constC'}
        Co = prop.C_J.*ones(ntime,1);
end
To = beta(:,:,1);
out = struct();

end



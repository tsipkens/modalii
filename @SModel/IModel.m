
% IMODEL  Evaluates the inverse spectroscopic model (using some form of pyrometry).
%  
%  T = SModel.IModel(PROP, J) uses the properties in PROP and incandescence
%  in J to compute temperatures for the given spectroscopic model. 
%  
%  NOTE: Uses SModel.opts.pyrometry to switch between evaluation methods.
%  NOTE: Additional outputs are available. 
%  
%  AUTHOR: Timothy Sipkens

function [To, Ti, Co, s_T, out] = IModel(smodel, prop, J)

nshots = length(J(1,:,1)); % number of shots
ntime = length(J(:,1,1)); % number of times

%-- Choose type of pyrometry ---------------------------------------------%
%   Note: Chooses between the different methods of pyrometry and output
%   parameter using the opts.pyrometry in smodel.
switch smodel.opts.pyrometry
    
    %-- Two-color pyrometry --%
    % Simple/fast two-color pyrometry, default if two wavelengths.
    % Does not currently output scaling factor. 
    case {'2color', 'ratio'}
        l = smodel.l;
        Emr = prop.Emr(l(1), l(2), prop.dp0); % two-colour pyrometry
        To = (0.0143877696*(1/(l(2)*1e-9)-1/(l(1)*1e-9)))./...
            log(J(:,:,1)./J(:,:,2).*(((l(1)/l(2))^6)/Emr)); % phi=0.01438
        To = real(To);
        s_T = [];
        out = [];
        Co = [];
        
    case {'2color-scalingfactor'}
        l = smodel.l;
        Emr = prop.Emr(l(1),l(2), prop.dp0); % two-colour pyrometry
        To = (0.0143877696*(1/(l(2)*1e-9)-1/(l(1)*1e-9)))./...
            log(J(:,:,1)./J(:,:,2).*(((l(1)/l(2))^6)/Emr)); % phi=0.01438
        To = real(To);
        Co = bsxfun(@times,J,1./smodel.FModel(prop, To, prop.Em));
        Co = Co(:,:,1);
        s_T = [];
        out = [];
        
    case {'2color-constT'}
        l = smodel.l;
        Emr = prop.Emr(l(1),l(2), prop.dp0); % two-colour pyrometry
        To = (0.0143877696*(1/(l(2)*1e-9)-1/(l(1)*1e-9)))./...
            log(J(:,:,1)./J(:,:,2).*(((l(1)/l(2))^6)/Emr)); % phi=0.01438
        To = real(To);
        Co = J./smodel.FModel(1730.*ones(size(To)),prop.Em);
        Co = Co(:,1,1);
        s_T = [];
        out = [];
        
    case {'2color-advanced'}  % calculate temperature, pre-averaged data
        data1 = mean(J(:,:,1),2);
        data2 = mean(J(:,:,2),2);
        [To,Co] = smodel.calcRatioPyrometry(data1,data2);
        
        nn = 1000; % used for sampling methods
        s1 = std(J(:,:,1),[],2)./sqrt(nshots);
        s2 = std(J(:,:,2),[],2)./sqrt(nshots);
        datas1 = (mvnrnd(data1,s1'.^2,nn))';
        datas2 = (mvnrnd(data2,s2'.^2,nn))';
        % s_C = 0.3;
        % Cs = normrnd(1,s_C,[1,nn]);
        % datas1 = bsxfun(@times,datas1,Cs);
        % datas2 = bsxfun(@times,datas2,Cs);
        [~,~,s_T,out] = smodel.calcRatioPyrometry(datas1,datas2);
        
        out.resid = zeros(ntime,1);
        
    case {'constC'} % hold C at some constant value (not working)
        if ~isfield(opts, 'C')
            error('Error occurred in pyrometry: C was not specified.');
            return
        end
        
        To = 1;
        s_T = [];
        out = [];
        Co = [];
        
    % Spectral fitting / inference **********************
    otherwise
        switch smodel.opts.multicolor
            case {'constC-mass','priorC-smooth'} % simultaneous inference
                [To,Co,s_T,out] = smodel.calcSpectralFit_all(J);
            otherwise % sequential inference
                [To,Co,s_T,out] = smodel.calcSpectralFit(J);
        end
end

Ti = nanmean(To(1,:)); % average temperatures

end



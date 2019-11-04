
function prop = N2(prop,opts)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('prop','var'); prop = struct(); end

if ~exist('opts','var'); opts = struct(); end
if ~isfield(opts,'propmdel'); opts.propmodel = 'default'; end
%-------------------------------------------------------------------------%


prop.mg = 4.65186e-26;
R = 8.314472; % ideal gas constant

switch opts.propmodel
    case 'Michelsen'
        cpg = @(T) R.*(...
            3.498./(T.^2).*exp(1./T)./((exp(1./T)-1).^2)+...
            0.98378.*(3353.3./T).^2.*exp(3353.3./T)./((exp(3353.3./T)-1).^2)+...
            T./38811);
        prop.gamma1 = @(T)cpg(T)./(cpg(T)-R);
    case 'Liu'
        prop.gamma1 = @(T)1.4221-1.8636e-4.*T+8.0784e-8.*T.^2-...
            1.6425e-11.*T.^3+1.2750e-15.*T.^4;
    case 'Kock'
        cpg = @(T)28.58+3.77e-3.*T-5e4./(T.^2);
        prop.gamma1 = @(T)cpg(T)./(cpg(T)-R);
    case 'Charwath'
        prop.gamma1 = @(T)1.3;
    case 'Will'
        prop.gamma1 = @(T)1.3009;
    otherwise % e.g. Sipkens, default
        prop.gamma1 = @(T)7/5; % = 1.4
end
prop.gamma2 = @(T)(prop.gamma1(T)+1)/(prop.gamma1(T)-1);

end


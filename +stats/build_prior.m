
% BUILD_PRIOR  Build a prior function.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-21

function pr = build_prior(x0, Lx, type)

if ~exist('type', 'var'); type = []; end
if isempty(type); type = 'default'; end

if ~exist('Lx', 'var'); Lx = []; end
if isempty(Lx); Lx = speye(length(x0)); end
if any(size(Lx) == 1); Lx = diag(1 ./ Lx); end  % standard deviations supplied


switch type
    case {'default','none','uninformed'} % do not use a prior
        pr = @(x) 1;
        
    case {'independent','vector-xpr','normal','vector'} % only apply priors on QoIs
        pr = @(x) Lx * (x - x0);
        
    case {'scalar'}
        pr = @(x) norm(Lx * (x - x0));
        
    otherwise
        disp('ERROR: Invalid prior scheme!');
        return;
end


end


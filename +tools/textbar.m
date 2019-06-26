%=== textbar.m ===========================================================%
% Function: Print out a text-based progress bar
% Author:   Samuel Grauer, 2017-11-16

function textbar(pct)
%-------------------------------------------------------------------------%
% Input:
%   pct     Progress                                 	[%]
%-------------------------------------------------------------------------%


%--- Initialization ------------------------------------------------------%
% Parameters
n_dot = 50;
n_str = 75;

% Parse input
if ~exist('pct','var'), pct = 0; end
if isempty(pct) || pct < 0, pct = 0;
elseif pct > 1, pct = 1; end
%-------------------------------------------------------------------------%


%--- Print progress ------------------------------------------------------%
if pct == 0
    str_out = ['[',repmat(' ',[1 n_dot]),'] 0.0%%'];
else
    fprintf(repmat(char(8),[1 n_str]));
    nc = ceil(pct*n_dot);
    
    str_b01 = repmat('-',[1 nc]);
    str_b02 = repmat(' ',[1 n_dot-nc]);
    str_out = ['[',str_b01,str_b02,'] ',num2str(100*pct,'%.0f'),'%%'];
end

fprintf([str_out,repmat(' ',[1 n_str-length(str_out)])]);
fprintf('\n');
%-------------------------------------------------------------------------%
end
%=== End =================================================================%
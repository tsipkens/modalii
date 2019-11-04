
% GET Get properties corresponding to provided strings
% Author: Timothy Sipkens, 2019-11-03
%=========================================================================%

function prop = get(strs,opts)

%-- Parse inputs ---------------------------------------------------------%
if ~iscell(strs); strs = {strs}; end

if ~exist('opts','var'); opts = struct(); end
if isempty(opts); opts = struct(); end
%-------------------------------------------------------------------------%


prop = struct(); % initialize prop
for ss=1:length(strs)
    prop = eval(['Prop.',strs{ss},'(prop,opts)']);
end

end


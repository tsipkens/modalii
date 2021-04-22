
% MLE2table  Convert MLE and fields to a table.
%  Useful for displaying results.
% 
%  AUTHOR: Timothy Sipkens, 2021-04-21

function t = mle2table(MLE, fields, varargin)

if size(MLE, 1) == 1
    MLE = MLE';
end

t = table(MLE, 'RowNames', fields);

if ~exist('SX', 'var'); SX = {}; end
for ii=1:2:length(varargin)
    t.(varargin{ii}) = varargin{ii + 1};
end

disp(' ');
disp(t);
disp(' ');

if nargout==0; clear t; end

end

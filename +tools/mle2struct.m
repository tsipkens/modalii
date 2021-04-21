
% MLE2STRUCT  Convert MLE and fields to a structure.
% 
%  AUTHOR: Timothy Sipkens, 2021-04-21

function s = mle2struct(mle, fields)

s = struct();
for ii=1:length(fields)
    s.(fields{ii}) = mle(ii);
end

end

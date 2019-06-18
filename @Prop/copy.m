function [copied] = copy(prop)
% creates of copy of the current object

copied = Prop; % initiate new signal

% copy all non-hidden properties
p = properties(prop);
for ii = 1:length(p)
    copied.(p{ii}) = prop.(p{ii});
end

end


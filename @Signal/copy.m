function [copied] = copy(signal,vec)
% creates of copy of the current object

copied = Signal(signal.l); % initiate new signal

% copy all non-hidden properties
p = properties(signal);
for ii = 1:length(p)
    copied.(p{ii}) = signal.(p{ii});
end

if nargin>1; % use vec to reduce number of wavelengths
    copied.data = copied.data(:,:,vec);
    copied.l = copied.l(vec);
end

end


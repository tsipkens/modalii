
% PARSE_SIGNAL Parses the typical signals structured array into its parts.
% Author: Timothy Sipkens, 2020-03-27
%=========================================================================%

function [s,t,l] = parse_signal(signal)

if ~isa(signal,'double')
    s = signal.data;
    t = signal.t;
    l = signal.l;
    
else % if is a double, just return the signal
    s = signal;
    t = [];
    l = [];
end

end


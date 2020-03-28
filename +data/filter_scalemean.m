
% Scale the data in the signal by the mean of all of the data.
% Promotes stability in inference of temperature.

function [s] = filter_scalemean(s)

s = s./mean(s(:));


end


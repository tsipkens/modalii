function [] = filter_scalemean(signal)
% Scale the data in the signal by the mean of all of the data.
% Promotes stability in inference of temperature.

signal.data = signal.data./mean(signal.data(:));


end


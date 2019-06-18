function [] = filter_noneg(signal)

signal.data(signal.data<0)=NaN;

end


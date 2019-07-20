function [signal] = loadFiles_ih(fname,l)

 % APB 2016, Ar-17, CO2-33, 49-N2, 43-He
 % APB 2015, Ar-14; Ne-71; N2O-67; N2-55; He-48; CO2-42; CO-32?
filename{1} = ['..\Data\',fname,...
    ' Ch1 Absolute Intensity.csv'];
filename{2} = ['..\Data\',fname,...
    ' Ch2 Absolute Intensity.csv'];

% Signal ******************************************************************
signal = Signal(l);
opts.t0 = 13; % first time in raw data to use (Sina data: 1, Mo: 13; Ag: 17)
opts.ta = 153; % artificial start time (normally 0)
opts.te = 603; % end time (Fe: 375; Mo: 603; Ag: 52)
opts.tt = 42;

disp('Loading signal...');
for ii=1:length(filename)
    signal.J_raw(:,:,ii) = csvread(filename{ii})';
end

signal.t_raw = signal.J_raw(:,1,1).*1e9; % copy times froms first column
signal.t = signal.t_raw; % make filter copy
signal.J_raw = signal.J_raw(:,2:end,:); % remove time column
signal.J = signal.J_raw;
signal.data = signal.J_raw;

signal.filter_tpos; % remove negative times

signal.filter_mismatch; % remove signals whose peaks are out of sync

signal.filter_trimnoise(opts); % trim noise at a certain threshold

signal.filter_ave(20); % Average every 20 signals
% signal.filter_ave_rand(20); % Average every 20 signals (randomly), not used
signal.outlier; % Thomson-tau

disp('Signal loaded.');

end


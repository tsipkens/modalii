function [] = plot(signal,varargin)
%Options: 'shots', 'subset', 'type={time,wavelength}', 'plotspec'

shots = 1:length(signal.data(1,:,1)); % by default, plot all wavelengths
type = 'time'; % by default, plot with time
subset = [];
plotspec = '.';

ii=1;
while ii<=(nargin-1);
    if strcmp(varargin{ii},'shots')
        shots = varargin{ii+1};
        ii = ii+2; % Skip an input
    elseif strcmp(varargin{ii},'type')
        type = varargin{ii+1};
        ii = ii+2; % Skip an input
    elseif strcmp(varargin{ii},'subset')
        subset = varargin{ii+1};
        ii = ii+2; % Skip an input
    elseif strcmp(varargin{ii},'plotspec')
        plotspec = varargin{ii+1};
        ii = ii+2; % Skip an input
    else
        disp(['Input ''',varargin{ii},''' is not an option.']);
        ii = ii+2; % Skip an input
    end
end

clf;
hold on;
switch type
    case 'time'
        for ii=shots;
            if isempty(subset); subset=1:length(signal.l); end;
            plot(signal.t,squeeze(signal.data(:,ii,subset)),plotspec,'markers',4);
        end
    case 'wavelength'
        for ii=shots;
            if isempty(subset); subset=1:length(signal.t); end;
            plot(signal.l,squeeze(signal.data(subset,ii,:))',plotspec,'markers',4);
        end 
end
hold off;

end

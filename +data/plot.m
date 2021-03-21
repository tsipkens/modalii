
function [] = plot(s,t_l,type,varargin)
%Options: 'shots', 'subset', 'type={time,wavelength}', 'plotspec'

if ~exist('type','var'); type = []; end
if isempty(type); type = 'time'; end % by default, plot with time


shots = 1:length(s(1,:,1)); % by default, plot all wavelengths
subset = [];
plotspec = '.';

ii=1;
while ii<=(nargin-1)
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
        for ii=shots
            if isempty(subset); subset=1:length(l); end
            plot(t_l,squeeze(s(:,ii,subset)),plotspec,'markers',4);
        end
    case 'wavelength'
        for ii=shots
            if isempty(subset); subset=1:length(t); end
            plot(t_l,squeeze(s(subset,ii,:))',plotspec,'markers',4);
        end 
end
hold off;

end

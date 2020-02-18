classdef Signal < handle
    
    properties
        data = []; % data to be used in inference
        
        t_raw = []; % stores raw times
        t = []; % stores time
        J_raw = []; % stores raw signal
        J = []; % stores processed signal *remove?
        
        l = []; % wavelength
        t0 = []; % starting time (processing var.)
        tp = []; % peak time (processing var.) *update to be a processing var.
        cut = []; % range of times used (processing var.)
        
        F0 = []; % signal fluence
        gas = []; % gas for signal
        matl = []; % matl for signal
        type = []; % origin of data
        
        opts struct = struct(...
            'range','default',...
            'isave',false); % options
    end
    
    methods
        function signal = Signal(l,varargin)
            signal.l = l;
            
            %-- Handle additional properties -----------------------------%
            ii = 1;
            while ii<=(nargin-1)
                if isprop(signal,varargin{ii})
                    signal.(varargin{ii}) = varargin{ii+1};
                    if strcmp(varargin{ii},'J_raw')
                        signal.J = signal.J_raw;
                    end
                    ii = ii+2;
                else
                    aa = fieldnames(varargin{ii});
                    bb = varargin{ii};
                    for jj = 1:length(aa)
                        if isfield(signal.opts,aa{jj})
                            signal.opts.(aa{jj}) = bb.(aa{jj});
                        end
                    end
                    ii = ii+1;
                end
            end
        end
        
        %-- Filters ------------------------------------------------------%
        [] = filter_ave(signal,nn); % Find mean of across given shots.
        [] = filter_ave_rand(signal,n); % Find mean of across given shots, randomize signals avereged.
        [] = filter_mismatch(signal,t0); % Remove signals with mismatched peaks.
        [] = filter_noneg(signal); % Replace negative values with NaN.
        [] = filter_norm(signal); % Normalize the signals by peak.
        [] = filter_null(signal); % Remove null entries. 
        [] = filter_switchneg(signal); % Invert signal and remove residual.
        [] = filter_sync(signal,varargin); % Sync signals according to peak.
        [] = filter_tpos(signal); % Remove non-increasing times.
        [] = filter_trimnoise(signal,t0); % Trim signal based on noise/peak. (could be updated)
        [] = filter_trimpeak(signal); % trim before the peak
        [] = filter_combinewave(signal,rr); % combine the measurement at several wavelengths to reduce noise
        [] = filter_signalmean(signal); % scale the data in the signal by the mean of all of the data
        [] = outlier(signal); % perform Thomson-Tau outlier removal on the data
        
        %-- Peak temperature calculations --------------------------------%
        [Tp,tp,F0] = get_peak_temp(signal,prop); % Calculate the peak temperature and time
        
        % Plotting functions ---------------------------------------------%
        [] = plot(signal,varargin); % Plot signals specified in num (or all, if num is excluded)
        [] = plot_image(signal,ind); % plot incandescence as an image
        
        
        [copied] = copy(signal,vec); % make a copy of the specified signal
    end
    
    methods(Static)
        [J] = blackbody(T,l,opts,varargin);
    end
    
end


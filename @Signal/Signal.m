
classdef Signal
    
    properties
        data = []; % data to be used in inference
        
        t = []; % stores time
        J = []; % stores processed signal *remove?
        J_raw = [];
        
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
    end
    
end


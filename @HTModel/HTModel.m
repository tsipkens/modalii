
classdef HTModel
% HTMODEL Class definition for the TiRe-LII heat transfer model.
    
    properties
        prop = []; % material and experimental properties
        t = []; % time vector
        
        x = []; % cell containing QoI variable names
        
        dTdt = []; % function handle, rate of temperature change
        dmdt = []; % function handle, rate of mass change
        dXdt = []; % function handle, rate of annealed fraction change
        
        opts@struct = struct(...
            'evap','default',... % evaporation model
            'cond','default',... % conduction model
            'rad','default',... % radiation model
            'abs','default',... % absorption model
            'ann','default',... % annealing model
            'laserprofile','default',... % temporal shape of laser profile
            'polydispersity',0,... % whether to incorporate polydispersity
            'deMethod','default'... % ODE sovler method
            );
    end
    
    methods
        %-- Constructor method -------------------------------------------%
        function htmodel = HTModel(prop,x,t,varargin)
            htmodel.prop = prop; % material and experimental properties
            htmodel.x = x; % cell containing QoI variable names
            htmodel.t = t; % time vector
            
            %-- Handle additional options --------------------------------%
            ii = 1;
            while ii<=length(varargin) % incorporate list of properties
                if isprop(htmodel,varargin{ii}) % manually set property
                    htmodel.(varargin{ii}) = varargin{ii+1};
                    ii = ii+2; % skip an input
                    
                else  % incorporate single opts struct
                    aa = fieldnames(varargin{ii});
                    bb = varargin{ii};
                    for jj = 1:length(aa)
                        if isfield(htmodel.opts,aa{jj})
                            htmodel.opts.(aa{jj}) = bb.(aa{jj});
                        end
                    end
                    ii = ii+1;
                    
                end
            end
            %-------------------------------------------------------------%
            
            htmodel = htmodel.gov_eqn; % generate governing equation
        end
        %-----------------------------------------------------------------%
        
        
        %-- Heat transfer solvers ----------------------------------------%
        [Tout] = evaluate(htmodel,x); % evaluates selected model at x, outputting temperature
        [Tout,dpo,mpo,Xo] = de_solve(htmodel,prop,dp) % solves ode at a specified particle size
        [dTdt,dmdt] = gov_eqn(htmodel) % determines governing equation and stores in dTdt
        %-----------------------------------------------------------------%
        
        
        %-- Heat transfer submodels --------------------------------------%
        [q] = q_cond(htmodel,T,dp); % Evaluates conduction at the specified parameters
        [q,J,hv,pv] = q_evap(htmodel,T,dp); % Evaluates evaporation at the specified parameters
        [J] = Jevap(htmodel,T,dp); % Evaluates rate of mass loss at the specified parameters
        [q,rad] = q_rad(htmodel,T,dp); % Evaluates radiation at the specified parameters
        [q,Cabs] = q_abs(htmodel,t,dp); % Evaluates the absorbed laser energy at specified parameters
        [q,dXdt] = q_ann_Mich(htmodel,T,dp,X); % Evaluates Michelsen's annealing model
        [q,dXdt] = q_ann_Sip(htmodel,T,dp,X); % Evaluates in house annealing model
        [dXdt] = dXdt_fun(htmodel,q_ann,T,dp,X); % Evaluates dXdt, second output from q_ann
        %-----------------------------------------------------------------%
        
        
        %-- Plotting functions -------------------------------------------%
        [q,h] = sig_q(htmodel,T,dp,opts_rad); % Evaluates the different modes at the specified parameters
        %-----------------------------------------------------------------%
        
    end
    
end


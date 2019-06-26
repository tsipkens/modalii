
% HTMODEL   Class definition for the TiRe-LII heat transfer model.
% Author:   Timothy Sipkens, 2015
%=========================================================================%

classdef HTModel
    
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
        % HTMODEL Constructor method for heat transfer model.
        %-----------------------------------------------------------------%
        % Inputs:
        %   prop        Instance of the properties class containing
        %               material properties
        %   x           Cell of strings defining fields of prop to be
        %                   considered as quantities of interest (Qoi)
        %                   (e.g. 'alpha','dp')
        %   t           Vector of times at which heat transfer model is
        %                   evaluated
        %   varargin    Options as a struct or a series of name-value pairs
        %                   (Optional)
        %-----------------------------------------------------------------%
        
            htmodel.prop = prop; % copy material and experimental properties
            htmodel.x = x; % copy cell containing QoI variable names
            htmodel.t = t; % time vector
            
            htmodel = tools.parse_varargin(htmodel,varargin{:});
                % handle additional options (see function in tools package)
            
            htmodel = htmodel.gov_eqn; % generate governing equation
        end
        
        
        %-- Heat transfer evaluation functions ---------------------------%
        [Tout] = evaluate(htmodel,x); % evaluates selected model at x, outputting temperature
        [Tout,dpo,mpo,Xo] = de_solve(htmodel,prop,dp) % solves ode at a specified particle size
        [dTdt,dmdt] = gov_eqn(htmodel) % determines governing equation and stores in dTdt
        %-----------------------------------------------------------------%
        
        
        %-- Heat transfer submodels --------------------------------------%
        [q] = q_cond(htmodel,T,dp); % evaluates conduction at the specified parameters
        [q,J,hv,pv] = q_evap(htmodel,T,dp); % evaluates evaporation at the specified parameters
        [J] = Jevap(htmodel,T,dp); % evaluates rate of mass loss at the specified parameters
        [q,rad] = q_rad(htmodel,T,dp); % evaluates radiation at the specified parameters
        [q,Cabs] = q_abs(htmodel,t,dp); % evaluates the absorbed laser energy at specified parameters
        [q,dXdt] = q_ann_Mich(htmodel,T,dp,X); % evaluates Michelsen's annealing model
        [q,dXdt] = q_ann_Sip(htmodel,T,dp,X); % evaluates in house annealing model
        [dXdt] = dXdt_fun(htmodel,q_ann,T,dp,X); % evaluates dXdt, second output from q_ann
        %-----------------------------------------------------------------%
        
        
        %-- Plotting functions -------------------------------------------%
        [q,h] = sig_q(htmodel,T,dp,opts_rad); % evaluates the different modes at the specified parameters
        %-----------------------------------------------------------------%
        
    end
    
end


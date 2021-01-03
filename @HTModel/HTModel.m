
% HTMODEL  Class definition for the TiRe-LII heat transfer model.
% AUTHOR: Timothy Sipkens, 2015
%=========================================================================%

classdef HTModel
    
    properties
        prop = [];  % default material and experimental properties
        t = [];  % time vector
        
        x = [];  % cell containing QoI variable names
        
        dTdt = [];  % function handle, rate of temperature change
        dmdt = [];  % function handle, rate of mass change
        dXdt = [];  % function handle, rate of annealed fraction change
        
        % Heat transfer model options.
        opts struct = struct(...
            'cond', 'free-molecular',...  % conduction model
            'evap', 'free-molecular',...  % evaporation model
            'rad', 'none',...  % radiation model
            'abs', 'none',...  % absorption model
            'ann', 'none',...  % annealing model
            'laserprofile', 'default',...  % temporal shape of laser profile
            'polydispersity', 0,...  % whether to incorporate polydispersity
            'deMethod', 'default'...  % ODE sovler method
            );
    end
    
    
    methods
        %== HTMODEL ======================================================%
        %   Constructor method for heat transfer model.
        % 
        % INPUTS:
        %   prop        Instance of the properties class containing
        %               material properties
        %   x           Cell of strings defining fields of prop to be
        %                   considered as quantities of interest (Qoi)
        %                   (e.g. 'alpha','dp')
        %   t           Vector of times at which heat transfer model is
        %                   evaluated
        %   varargin    Options as a struct or a series of name-value pairs (Optional)
        %-----------------------------------------------------------------%
        function htmodel = HTModel(prop, x, t, varargin)
            
            % If no fields specified, only use particle diameter.
            if isempty(x); x = {'dp0'}; end
            
            htmodel.prop = prop; % copy material and experimental properties
            htmodel.x = x; % copy cell containing QoI variable names
            htmodel.t = t; % time vector
            
            % Handle additional options (see function in tools package).
            htmodel = tools.parse_varargin(htmodel, varargin{:});
            
            % Print HTModel properties to console.
            tools.textheader('New heat transfer model');
            disp(['  Conduction:  ', htmodel.opts.cond]);
            disp(['  Evaporation: ', htmodel.opts.evap]);
            disp(['  Absorption:  ', htmodel.opts.abs]);
            disp(['  Radiation:   ', htmodel.opts.rad]);
            disp(['  Annealing:   ', htmodel.opts.ann]);
            tools.textheader();
            
        end
        
        
        %-- Heat transfer evaluation functions ---------------------------%
        [Tout] = evaluate(htmodel, x); % evaluates selected model at x, outputting temperature
        [Tout, dpo, mpo, Xo] = de_solve(htmodel, prop, dp); % solves ode at a specified particle size
        [htmodel, dTdt, dmdt] = de_build(htmodel); % determines governing equation and stores in dTdt
        %-----------------------------------------------------------------%
        
        
        %-- Heat transfer submodels --------------------------------------%
        [q] = q_cond(htmodel, prop, T, dp); % evaluates conduction at the specified parameters
        [q, J, hv, pv] = q_evap(htmodel, prop, T, dp); % evaluates evaporation at the specified parameters
        [J] = Jevap(htmodel, prop, T, dp); % evaluates rate of mass loss at the specified parameters
        [q, rad] = q_rad(htmodel, prop, T, dp); % evaluates radiation at the specified parameters
        [q, Cabs] = q_abs(htmodel, prop, T, dp); % evaluates the absorbed laser energy at specified parameters
        [q, dXdt] = q_ann_Mich(htmodel, prop, T, dp,X); % evaluates Michelsen's annealing model
        [q, dXdt] = q_ann_Sip(htmodel, prop, T, dp,X); % evaluates in house annealing model
        [dXdt] = dXdt_fun(htmodel,q_ann,T,dp,X); % evaluates dXdt, second output from q_ann
        %-----------------------------------------------------------------%
        
        
        %-- Plotting functions -------------------------------------------%
        [q, h] = sig_q(htmodel,T,dp,opts_rad); % evaluates the different modes at the specified parameters
        %-----------------------------------------------------------------%
        
    end
    
end


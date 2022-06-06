
% HTMODEL  Class definition for the TiRe-LII heat transfer model.
%  The htmodel.evaluate(...) function is then to be used to evaluate the 
%  heat transfer model for the fields specified by x. 
%  
%  htmodel = HTModel(PROP, [], t) creates a heat transfer model using the
%  default options, the material and experimental properties in PROP, 
%  and at the times specified by vector t. The model is to be evaluated by 
%  providing a primary particle diameter (i.e., this is equivalent to 
%  X = 'dp0') using htmodel.evaluate(dp0). 
%  
%  htmodel = HTModel(PROP, X, t) evaluates the heat transfer model, as
%  above, but allows for evaluation at the properties named in X. For
%  example, X = {'dp0', 'alpha'} will create a model that is a funciton of
%  primary particle size and thermal accommodation coefficient (e.g.,
%  htmodel.evaluate([30, 0.5]) will then evaluate the heat transfer model
%  at dp0 = 30 and alpha = 0.5. Most fields in PROP can be specified in X.
%  
%  htmodel = HTModel(PROP, X, t, OPTS) add an options structure that
%  controls the nature of the heat transfer model. For example, 
%  OPTS.abs = 'include' will add an laser absorption mode to the model. For
%  the list of fields, see the class properties in the HTModel.m file.
%  Options can also be name-value pairs. 
%  
%  AUTHOR: Timothy Sipkens, 2015

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
        %   See class header for use.
        function htmodel = HTModel(prop, x, t, varargin)
            
            % If no fields specified, only use particle diameter.
            if ~exist('x', 'var'); x = []; end
            if isempty(x); x = {'dp0'}; end
            if ~iscell(x); x = {x}; end
            
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
        [T, mpo, Xo] = evaluate(htmodel, x);
        [Tout, dpo, mpo, Xo] = de_solve(htmodel, prop, dp); % solves ode at a specified particle size
        [htmodel, dTdt, dmdt] = de_build(htmodel); % determines governing equation and stores in dTdt
        %-----------------------------------------------------------------%
        
        
        %-- Heat transfer submodels --------------------------------------%
        [q, Kn] = q_cond(htmodel, prop, T, dp, opts_cond); % evaluates conduction at the specified parameters
        [q, J, hv, pv] = q_evap(htmodel, prop, T, dp); % evaluates evaporation at the specified parameters
        [J] = Jevap(htmodel, prop, T, dp); % evaluates rate of mass loss at the specified parameters
        [q, rad] = q_rad(htmodel, prop, T, dp); % evaluates radiation at the specified parameters
        [q, Cabs] = q_abs(htmodel, prop, T, dp); % evaluates the absorbed laser energy at specified parameters
        [q, dXdt] = q_ann_Mich(htmodel, prop, T, dp,X); % evaluates Michelsen's annealing model
        [q, dXdt] = q_ann_Sip(htmodel, prop, T, dp,X); % evaluates in house annealing model
        [dXdt] = dXdt_fun(htmodel,q_ann,prop,T,dp,X); % evaluates dXdt, second output from q_ann
        %-----------------------------------------------------------------%
        
        
        %-- Plotting functions -------------------------------------------%
        [q, h] = sig_q(htmodel,T,dp,opts_rad); % evaluates the different modes at the specified parameters
        %-----------------------------------------------------------------%
        
    end
    
end


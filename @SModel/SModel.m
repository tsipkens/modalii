
% SMODEL  Class for the spectroscopic model and pyrometry calculations.
%  
%  smodel = SMODEL(PROP, X, t, L) generates a spectroscopic model for the
%  parameters given in the cell, X, and for the material properties in
%  PROP, the times in t, and the wavelengths in L. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2017

classdef SModel
    
    properties
        prop = [];
        l = []; % wavelength
        t = []; % time
        T = []; % temperature, as function handle
        J = []; % incandescence
        htmodel = []; % embedded heat transfer model
        
        x = []; % variable names, quantities of interest (QoIs)
        
        data_sc = []; % used to scale Planck's to stabalize inference algorithms
        
        opts struct = struct( ...
            'multicolor', 'default', ... % indicates which multicolor sovler to use
            'pyrometry', 'ratio' ... % indicates how to handle pyrometry
            );
    end
    
    methods
        %-- Constructor method -------------------------------------------%
        function smodel = SModel(prop, x, t, l, varargin)
            smodel.prop = prop;
            smodel.x = x;
            smodel.t = t;
            smodel.l = l;
            
            % Handle additional options (see function in tools package).
            smodel = tools.parse_varargin(smodel, varargin{:});
            
            T_sc = 3000; % temperature used for scaling/stability, [K]
            smodel.data_sc = smodel.blackbody(T_sc,1064).*...
                prop.Em(1064,30)/(1064e-9); % scale Planck's law about T_sc for stability
                    % use dp = 30 nm for data scaling (used for stability)
                    
            tools.textheader('New spectroscopic model');
            disp(['  Pyrometry: ', smodel.opts.pyrometry]);
            tools.textheader();
        end
        %-----------------------------------------------------------------%
        
        %-- Modeling functions -------------------------------------------%
        [Jo] = FModel(smodel, prop, T, Em) % calculates J given T, Em is function handle
        [To,Ti,Co,s_T,s_C,r_T,resid,oth] = IModel(smodel, prop, J) % solve inverse model for temperature
        [To,Co,s_T,out] = calcSpectralFit(smodel,J) % spectral fitting with sequential inference
        [To,Co,s_T,out] = calcSpectralFit_all(smodel,J) % spectral fitting with simulatneous inference
        [To,Co,s_T,out] = calcRatioPyrometry(smodel,J1,J2,Emr) % ratio pyrometry evaluation, corr. calc. included
        
        %-- Evaluate functions, take and update x ------------------------%
        [To,Co] = evaluateI(smodel, x, J) % evaluate inverse model at x, given J in SModel
        [To,Jo] = evaluateIF(smodel, x) % evaluate inverse and forward model at x, given T in Smodel
        [Jo,mp] = evaluateF(smodel, x) % evaluate forward model at x
        
        %-- Plotting functions -------------------------------------------%
        [] = plotI(smodel,T,C,J,n) % Plots data vs fit of inverse procedure
        
    end
    
    methods(Static)
        [J] = blackbody(T,l);
        [b] = outlier(b); % remove outlier by Thomson-Tau procedure
    end
    
end



% SMODEL    Class for the spectroscopic model and pyrometry calculations.
% Author:   Timothy Sipkens, 2017
%=========================================================================%

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
        
        opts@struct = struct( ...
            'multicolor','default', ... % indicates which multicolor sovler to use
            'pyrometry','default' ... % indicates how to handle pyrometry
            );
    end
    
    methods
        %-- Constructor method -------------------------------------------%
        function smodel = SModel(prop,x,t,l,varargin)
            smodel.prop = prop;
            smodel.x = x;
            smodel.t = t;
            smodel.l = l;
            
            %-- Parse additional innputs ---------------------%
            ii = 1;
            while ii<=length(varargin)
                if isa(varargin{ii},'Signal') % derive paramters from signal
                    smodel.J = varargin{ii}.data;
                    ii = ii+1;
                    
                elseif isa(varargin{ii},'HTModel') % derive paramters from heat transfer model
                    smodel.htmodel = varargin{ii};
                    ii = ii+1;
                    
                elseif isprop(smodel,varargin{ii}) % manually set property
                    smodel.(varargin{ii}) = varargin{ii+1};
                    ii = ii+2; % skip an input
                    
                else  % incorporate opts variable
                    aa = fieldnames(varargin{ii});
                    bb = varargin{ii};
                    for jj = 1:length(aa)
                        if isfield(smodel.opts,aa{jj})
                            smodel.opts.(aa{jj}) = bb.(aa{jj});
                        end
                    end
                    ii = ii+1;
                    
                end
            end
            
            T_sc = 3000; % temperature used for scaling/stability, [K]
            smodel.data_sc = smodel.blackbody(T_sc,1064).*...
                prop.Em(1064,30)/(1064e-9); % scale Planck's law about T_sc for stability
                    % use dp = 30 nm for data scaling (used for stability)
        end
        %-----------------------------------------------------------------%
        
        %-- Modeling functions -------------------------------------------%
        [Jout] = FModel(smodel,T,Em) % calculates J given T, Em is function handle
        [Tout,Ti,Cout,s_T,s_C,r_T,resid,oth] = IModel(smodel,J) % solve inverse model for temperature
        [Tout,Cout,s_T,out] = calcSpectralFit(smodel,J) % spectral fitting with sequential inference
        [Tout,Cout,s_T,out] = calcSpectralFit_all(smodel,J) % spectral fitting with simulatneous inference
        [Tout,Cout,s_T,out] = calcRatioPyrometry(smodel,J1,J2) % ratio pyrometry evaluation, corr. calc. included
        
        %-- Evaluate functions, take and update x ------------------------%
        [Tout] = evaluateI(smodel,x) % evaluate inverse model at x, given J in SModel
        [Tout] = evaluateIF(smodel,x) % evaluate inverse and forward model at x, given T in Smodel
        [Jout,mp] = evaluateF(smodel,x) % evaluate forward model at x
        
        %-- Plotting functions -------------------------------------------%
        [] = plotI(smodel,T,C,J,n) % Plots data vs fit of inverse procedure
        
    end
    
    methods(Static)
        [J] = blackbody(T,l);
        [b] = outlier(b); % remove outlier by Thomson-Tau procedure
    end
    
end


classdef Stats < handle
    
    properties
        model = []; % model as a function handle, can be non-linear
        b = []; % data vector, can be function handle specified or vector
        
        %-- Data noise --%
        Gb = [] % noise information, covariance matrix
        Gb_i = []; % noise information, inverse covariance matrix
        Lb = []; % transpose of Cholesky factorization of inverse covariance matrix
        sb = []; % noise information, standard deviation of each pt.
        
        %-- QoIs --%
        Gx_pr = []; % prior covariance matrix
        Lx_ipr = []; % chol. fact. of inverse prior covariance matrix
        x_pr = []; % prior expected value
        sx_pr = []; % prior std. dev. of QoIs
        
        %-- Nuisance parameters --%
        Gt_pr = []; % prior covariance matrix
        Lt_ipr = []; % chol. fact. of inverse prior covariance matrix
        t_pr = []; % prior expected value
        st_pr = []; % prior std. dev. of nuisance parameters
        
        pr_fun = []; % prior function
        min_fun = []; % minimization function
        
        opts@struct = struct(...
            'variance','default',... % method for using/calculating covariance in data
            'minimize','default',... % indicates minimization scheme
            'bFun',1,... % boolean indicating if data is dependent on x
            'prior','default',... % type of prior, default is joint normal
            'display','default'... % type of display in minimize function
            ); % Flags
    end
    
    methods
        %== Constructor method ===========================================%
        function stats = Stats(model,b,varargin)
            
            stats = stats.setup_b(b); % handle provided data
            
            %-- Parse inputs ---------------------------------------------%
            ii = 1;
            while ii<=length(varargin)
                if isprop(stats,varargin{ii}) % manually set property
                    stats.(varargin{ii}) = varargin{ii+1};
                    ii = ii+2; % skip an input
                else  % incorporate opts variable
                    aa = fieldnames(varargin{ii});
                    bb = varargin{ii};
                    for jj = 1:length(aa)
                        if isfield(stats.opts,aa{jj})
                            stats.opts.(aa{jj}) = bb.(aa{jj});
                        end
                    end
                    ii = ii+1;
                end
            end
            
            %-- Assign remaining properties ------------------------------%
            stats.model = model;
            
        end
        %=================================================================%
        
        
        %== SETUP_B ======================================================%
        %   Converts b to an appropriate form used elsewhere in the class
        %   Author:  Timothy Sipkens, 2019-07-18
        function stats = setup_b(stats,b)
        
        if ~isa(b,'function_handle')
            nn1 = size(b,1);
            nn2 = size(b,2);

            switch stats.opts.variance
                case {'covariance','dependent'}
                    Gb_i = inv(nancov(real(b'))./(nn2-1));
                    Lb = chol(stats.Gb_i);
                    bout = nanmean(b,2);
                case {'default','variance','independent'}
                    sb = nanstd(b,0,2)./sqrt(nn2-1);
                    Gb_i = spdiags(1./sb.^2,0,nn1,nn1);
                    Lb = chol(Gb_i);
                    bout = nanmean(b,2);
                case 'specified'
                    bout = nanmean(b,2);
                case 'single-shot'
                    sb = (0.05*mean(b)).*...
                        ones(nn1,nn2);
                    Gb_i = spdiags(1./sb.^2,...
                        0,nn1,nn1);
                    Lb = sqrt(Gb_i);
                    bout = b;
            end
            
            stats.Gb_i = @(x) Gb_i;
            stats.Lb = @(x) Lb;
            stats.b = @(x) bout;
        else
            nn1 = @(x)size(b(x),1);
            nn2 = @(x)size(b(x),2);

            switch stats.opts.variance
                case 'covariance'
                    Gb_i = @(x) inv(nancov(real(b(x)'))./...
                        (nn2(x)-1));
                    Lb = @(x) chol(Gb_i);
                    bout = @(x)nanmean(b(x),2);
                case {'default','variance','independent'}
                    sb = @(x) nanstd(b(x),0,2)./sqrt(nn2(x)-1);
                    Gb_i = @(x) spdiags(1./sb(x).^2,0,nn1(x),nn1(x));
                    Lb = @(x) chol(Gb_i(x));
                    bout = @(x) nanmean(b(x),2);
                case 'specified'
                    bout = @(x) nanmean(b(x),2);
                    Lb = [];
                    sb = [];
                    Gb_i = [];
                case 'single-shot'
                    sb = @(x) (0.05*mean(b(x))).*...
                        ones(nn1(x),nn2(x));
                    Gb_i = @(x) spdiags(1./sb(x).^2,...
                        0,nn1(x),nn1(x));
                    Lb = @(x) sqrt(Gb_i(x));
                    bout = b;
            end
            
            stats.Gb_i = Gb_i;
            stats.Lb = Lb;
            stats.b = bout;
            if exist('sb','var'); stats.sb = sb; end
        end

        end
        %=================================================================%
        
        
        %== Oher functions ===============================================%
        [min_fun] = get_min_fun(stats) % generate appropriate minimization function
        [pr_fun,pr_fun_theta] = get_prior_fun(stats) % generate prior functions
        [mle,jcb] = minimize(stats,x0,opts_or) % find MLE or MAP
        [post,x,y] = eval_post(stats,X,Y) % evaluate the posterior on the specified grid (ignore nuisance parameters)
        
        %-- Used in credible interval estimation -------------------------%
        [rnd,Gx,R,s] = cred_mcmc(stats,n,mle,sx_pr,theta) % get credible intervals from sampling
        [Gx_po,R,s] = cred_linear(stats,jcb) % get credible intervals from Jacobian
        [smp,G,R,s] = cred_bootstrap(stats,n,mle) % get credible intervals from bootstrapping
        [jcb] = jacob_est(stats,mle,min_fun) % estimate the jacobian about a point
        
        %-- Plotting functions -------------------------------------------%
        [] = plot_mle(stats,mle,t) % plot b_mod and b_exp
        [X,Y,x,y] = gen_grid(stats,mle,s,n) % generate grid to plot on
    end
    
end


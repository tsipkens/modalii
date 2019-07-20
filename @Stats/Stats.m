classdef Stats < handle
    
    properties
        model = []; % model as a function handle, can be non-linear
        b = []; % data vector, can be function handle specified by flag
        
        % Data noise ***********************
        Gb = [] % noise information, covariance matrix
        Gb_i = []; % noise information, inverse covariance matrix
        Lb_i = []; % Cholesky factorization of inverse covariance matrix
        sb = []; % noise information, standard deviation of each pt.
        
        % QoIs ***************************
        Gx_pr = []; % prior covariance matrix
        Lx_ipr = []; % chol. fact. of inverse prior covariance matrix
        x_pr = []; % prior expected value
        sx_pr = []; % prior std. dev. of QoIs
        
        % Nuisance parameters *************
        Gt_pr = []; % prior covariance matrix
        Lt_ipr = []; % chol. fact. of inverse prior covariance matrix
        t_pr = []; % prior expected value
        st_pr = []; % prior std. dev. of nuisance parameters
        
        pr_fun = []; % prior function
        min_fun = []; % minimization function
        
        opts@struct = struct(...
            'variance','default',... % method for using/calculating covariance in data
            'minimize','default',... % indicates minimization scheme
            'bFun',0,... % boolean indicating if data is dependent on x
            'prior','default',... % type of prior, default is joint normal
            'display','default'... % type of display in minimize function
            ); % Flags
    end
    
    methods
        function stats = Stats(model,b,varargin)
            stats.model = model;
            
            % parse additional inputs
            ii = 1;
            while ii<=length(varargin);
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
            
            % analyze b input
            if isa(b,'function_handle')
                stats.opts.bFun = 1; % make b into a static function
                % add determination of length of data produced by function
                %   handle?
            end
            
            % if the covariance was already specified
            if ~isempty(stats.Lb_i)
                stats.b = nanmean(b,2);
            else
                [stats.b,stats.Lb_i,stats.sb,stats.Gb_i] = stats.getNoiseInfo(b);
            end
            
        end
        
        [bout,Lb_i,sb,Gb_i] = getNoiseInfo(stats,bin) % gets mean and variance/covariance of data for analysis
        
        [min_fun] = getMinFun(stats) % generate appropriate minimization function
        [pr_fun,pr_fun_theta] = getPriorFun(stats) % generate prior functions
        [mle,jcb] = minimize(stats,x0,opts_or) % find MLE or MAP
        [post,x,y] = evalPost(stats,X,Y) % evaluate the posterior on the specified grid (ignore nuisance parameters)
        
        % Used in credible interval estimation **********************
        [rnd,Gx,R,s] = credMCMC(stats,n,mle,sx_pr,theta) % get credible intervals from sampling
        [Gx_po,R,s] = credLinear(stats,jcb) % get credible intervals from Jacobian
        [smp,G,R,s] = credBootstrap(stats,n,mle) % get credible intervals from bootstrapping
        [jcb] = jcbEst(stats,mle,min_fun) % estimate the jacobian about a point
        
        % Plotting functions ****************************************
        [] = plotmle(stats,mle,t) % plot b_mod and b_exp
        [X,Y,x,y] = genGrid(stats,mle,s,n) % generate grid to plot on
    end
    
end


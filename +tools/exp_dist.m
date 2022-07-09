
% EXP_DIST  Regularization based on the exponential of the distance between elements/pixels.
%  
%  X = invert.exp_dist(A, B, LAMBDA, LPR) performs exponential distance
%  regularization of the basic problem AX = B, solving for X. The quantity
%  LAMBDA is a regularization parameter, and pre-computed LPR is the 
%  structure of the exponential distance prio matrix. 
%  See invert.exp_dist_lpr(...) for more information on LPR. 
%  
%  X = invert.exp_dist(A, B, LAMBDA, LD, VEC) replaces the pre-computed
%  matrix with a length scale, LD, and vector at which the size
%  distribution is to be reconstructed. LPR is then built in the function. 
%  
%  X = invert.exp_dist(A, B, LAMBDA, LPR, [], XI) adds an initial X for
%  iterative solvers. {LPR, []} can be replaced with {LD, VEC}. 
%  
%  X = invert.exp_dist(..., SOLVER) adds a string to specify the type of
%  solver to use. 
%  
%  ------------------------------------------------------------------------
%
%  OUTPUTS:
%   X           Regularized estimate
%   D           Inverse operator (x = D*[b;0])
%   LPR0        Cholesky factorization of prior covariance
%   GPO_INV     Inverse of the posterior covariance after inversion
%  
%  ------------------------------------------------------------------------
%  
%  Author: Timothy Sipkens, 2018-10-22

function [x,D,Lpr0,Gpo_inv] = exp_dist(A, b, lambda, ld, vec, xi, solver)

x_length = size(A,2);

%-- Parse inputs -----------------------------------------------%
if ~exist('xi','var'); xi = []; end % if no initial x is given
if ~exist('solver','var'); solver = []; end % if computation method not specified

Lpr0 = exp_dist_lpr(ld, vec); % use external function to evaluate prior covariance
%--------------------------------------------------------------%

% Scael Lpr by regularization parameter.    
Lpr = lambda .* Lpr0;


%-- Choose and execute solver --------------------------------------------%
[x,D] = lsq(...
    [A;Lpr], [b;sparse(x_length,1)], xi, solver);


%-- Uncertainty quantification -------------------------------------------%
if nargout>=4
    Gpo_inv = A'*A+Lpr'*Lpr;
end

end




% EXP_DIST_LPR  A helper function for 'exp_dist' to compute prior covariance.
%  
%  LPR = invert.exp_dist_lpr(LD, VEC) computes the prior covariance using a
%  length scale of LD and a vector of sizes on which the distribution is to
%  be reconstructed, VEC. 
%  
%  LPR = invert.exp_dist_lpr(..., BC) adds an input to specify boundary
%  conditions. Currently, BC = 0 pushes for 0s that boundary.
%  
%  ------------------------------------------------------------------------
%  
%  Author: Timothy Sipkens, 2019-12-11

function [Lpr, D, Gpr] = exp_dist_lpr(ld, vec, bc)

%-- Parse inputs ------%
if ~exist('bc', 'var'); bc = []; end
if isempty(bc); bc = 1; end


%-- Compute distances between elements -----------------------------------%
[vec_a, vec_b] = ndgrid(vec,vec);  % grid of distances

dr = log10(vec_a) - log10(vec_b);
D = abs(dr ./ ld); % distance, abs replaces sqrt(d.^2)


%-- Compute prior covariance matrix --------------------------------------%
Gpr = exp(-D);

Gpr_inv = pinv(Gpr);
[Lpr, ~] = chol(Gpr_inv);
clear Gpr_inv; % to save memory

% NOTE: Will cause incompatibility between Lpr and Gpr.
if bc == 0  % then push for 0 at boundary
    Lpr(1, :) = 0;
    Lpr(1, 1) = Lpr(2, 2);
    Lpr(end, :) = 0;
    Lpr(end, end) = Lpr(2, 2);
end

% Remove entries where relative distance to current 
% pixel is greater than some constant. 
Lpr(D > 1.75) = 0;
Lpr = sparse(Lpr);

end





% LSQ      Performs least-squares or equivalent optimization.
%          For use with inversion functions like tikhonov.m
%          and exp_dist.m
% Author:  Timothy Sipkens, 2019-07-17
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   xi      Initial guess for solver    (Optional, default is empty)
%   solver  Solver                      (Optional, default is interior-point)
%
% Outputs:
%   x       Regularized estimate
%   D       Inverse operator (x = D*[b;0])
%=========================================================================%

function [x,D] = lsq(A,b,xi,solver)


%-- Parse inputs ---------------------------------------------------------%
if ~exist('xi','var'); xi = []; end % use empty if not specified
if ~exist('solver','var'); solver = []; end

if isempty(solver); solver = 'interior-point'; end % if computation method not specified
%-------------------------------------------------------------------------%


%-- Perform least-squares ------------------------------------------------%
x_length = length(A(1,:));
x_lb = sparse(x_length,1); % enforce non-negativity
if ~isempty(xi)
    xi = max(xi,x_lb); % modify to accommodate bound
end

%== Choose a solver and evaluate =========================================%
switch solver
    case 'non-neg' % constrained, iterative linear least squares
        options = optimset('Display','off','TolX',eps*norm(A,1)*length(A));
        x = lsqnonneg(A,b,options);
        D = []; % not specified when using this method
        
    case 'interior-point' % constrained, iterative linear least squares
        options = optimoptions('lsqlin','Algorithm','interior-point','Display','none');
        x = lsqlin(A,b,...
            [],[],[],[],x_lb,[],xi,options);
        D = []; % not specified when using this method
        
    case 'trust-region-reflective'
        D = (A'*A)\A'; % invert combined matrices to get first guess
        x_lb = D*b; % Note: previously, [b;Lx*zeros(x_length,1)]
        
        options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
        x = lsqlin(A,b,...
            [],[],[],[],x_lb,[],xi,options);
        
    case 'interior-point-neg'
        options = optimoptions('lsqlin','Algorithm','interior-point','Display','none');
        x = lsqlin(A,b,...
            [],[],[],[],[],[],xi,options);
        D = []; % not specified when using this method
        
    case 'algebraic' % matrix multiplication least squares (not non-negative constrained)
        D = (A'*A)\A'; % invert combined matrices
        x = D*b;
        
    case 'algebraic-inv' % alternate algebraic least squares (less stable than previous option)
        D = inv(A'*A)*A'; % invert combined matrices
        x = D*b;
        
end
%-------------------------------------------------------------------------%


end


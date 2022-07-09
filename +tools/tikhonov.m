
% TIKHONOV  Performs inversion using various order Tikhonov regularization in 2D.
%
%  INPUTSs:
%   A           Model matrix
%   b           Data
%   order       Order of regularization
%   xi          Initial guess for solver
%                   (OPTIONAL, default is zeros)
%   f2          Flag of whether to use two-stage procedure of 
%               Huckle and Sedlacek.
%   solver      Solver (OPTIONAL, default is 'interior-point')
%
%  OUTPUTS:
%   x           Regularized estimate
%   D           Inverse operator (x = D*[b;0])
%   Lpr0        Tikhonov matrix structure
%   Gpo_inv     Inverse of posterior covariance
%  
%  AUTHOR: Timothy Sipkens, 2020-04-11

function [x,D,Lpr0,Gpo_inv] = tikhonov(A, b, lambda, order, bc, xi, solver)

n = size(A,2);  % length of x

%-- Parse inputs ---------------------------------------------------------%
if ~exist('order', 'var'); order = []; end
    % if order not specified, use default of tikhonov_lpr

if ~exist('xi', 'var'); xi = []; end % if initial guess is not specified
if ~exist('solver', 'var'); solver = []; end

if ~exist('f2', 'var'); f2 = []; end
if isempty(f2); f2 = []; end
%-------------------------------------------------------------------------%


%-- Get Tikhonov smoothing matrix ----------------------------------------%
Lpr0 = tikhonov_lpr(...
    order, n, bc);
Lpr = lambda .* Lpr0;


%-- Choose and execute solver --------------------------------------------%
pr_length = size(Lpr0,1);
[x,D] = lsq(...
    [A;Lpr], [b;sparse(pr_length,1)], ...
    xi, solver);


%-- Uncertainty quantification -------------------------------------------%
if nargout>=4
    Gpo_inv = A'*A+Lpr'*Lpr;
end

end



% TIKHONOV_LPR  Generates Tikhonov matrices/operators. 
%  
%  L = tikhonov_lpr(ORDER, X_LENGTH) generates the Tikhonov matrix, L, of
%  the order specified in ORDER and corresponding to a vector x of length
%  X_LENGTH, such that L*x = 0.
%  
%  L = tikhonov_lpr(ORDER, X_LENGTH, BC) add an input for the type of
%  boundary conidition to be applied at the edges.
%  
%  ------------------------------------------------------------------------
%
%  ORDER options:
%   0: 0th order Tikhonov promotes small solutions. 
%   1: 1st order Tikhonov minimizes the first derivative.
%   2: 2nd order Tikhonov minimized the second derivative.
%  
%  BC (boundary condition) options:
%   []: No explicit boundary conditions (uses that implied by order).
%    0: Zero BCs.
%    1: No slope BCs
%    2: No curvature BCs.
%  
%  ------------------------------------------------------------------------
% 
%  AUTHOR:   Timothy Sipkens, 2020-04-11

function L = tikhonov_lpr(order, x_length, bc)

% Parse type input.
if ~exist('bc', 'var'); bc = []; end
if isempty(bc); bc = order; end  % neglects boundary conditions

% Choose between order of Tikhonov operator to generate.
switch order
    case 0 % 0th order Tikhonov
    	L = speye(x_length);
        
    case 1 % 1st order Tikhonov
        L = -speye(x_length);
        L = spdiags(ones(x_length, 1), 1, L);
        
        if bc~=0
            L(end,:) = [];
        end
        
    case 2 % 2nd order Tikhonov
        L = -speye(x_length);
        L = spdiags(0.5 .* ones(x_length,2), [-1,1], L);
        
        if bc == 0  % zero
            L(1, 2) = 0;
            L(x_length, x_length - 1) = 0;
       elseif bc == 1  % no slope
            L(1, 2) = 1;
            L(x_length, x_length - 1) = 1;
        else
            L(1,:) = [];
            L(end,:) = [];
        end
        
    otherwise
        disp('The specified order of Tikhonov is not available.');
        disp(' ');
        return
end

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


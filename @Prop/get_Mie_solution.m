
function [Em_eff,Q_abs] = get_Mie_solution(n,k,lambda_vec,dp_vec)
% Get Mie solution from MatScat code. 

% lambda_vec = 442e-9;%(300:1:800).*1e-9; % range of wavelengths to be considered
% dp_vec = 20e-9;%(10:1:100).*1e-9 ; % range of nanosphere diameter

% [~,~,nAu] = LD(lambda_vec*1e-9,'Au','LD'); %Drude model for Liquid Si

if ~isa(n,'function_handle') % check if inputs are function handles
    n_C = n;
    n = @(l) n_C; % convert ot function handle
    k_C = k;
    k = @(l) k_C; % convert ot function handle
end

nm = 1.0;           % outer medium refractive index (real)
nang = 1800;        % number of far field angles to evaluate
conv = 1;   % convergence factor

Q_abs = zeros(length(dp_vec),length(lambda_vec));
Em_eff = Q_abs;

for jj = 1:length(dp_vec)
    dp = dp_vec(jj);
    for kk = 1:length(lambda_vec)
        
        lambda = lambda_vec(kk); 
        ns = n(lambda)+k(lambda)*1i;    % sphere refractive index (complex)

        rad = dp/2;           % sphere radius
        wavenumber = 2*pi/lambda*nm;    % wavenumber in medium n_m
        
        % Calculate cross sections *****************
        [~, C, ~] = MatScat_Mie_code.calcmie(rad, ns, nm, ...
            lambda, nang, 'ConvergenceFactor', conv);
       
        % Calculate cross sections and efficiencies ************
        Q = MatScat_Mie_code.util.getEfficiencies(C,rad(end),3);
        % Q_abs(jj,kk) = Q.abs/(pi*(dp/2)^2);
        Q_abs(jj,kk) = Q.abs;
        Em_eff(jj,kk) = Q.abs/(4*pi*dp/lambda);
    end
end

end



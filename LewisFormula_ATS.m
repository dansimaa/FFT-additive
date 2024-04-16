function prices = LewisFormula_ATS(moneyness, F0, B, t, alpha, ATS_params)

% The function computes the prices of European call options, given as input their log-moneyness,
% using the Lewis formula with the characteristic function of an  ATS process with power-law scaling 
% parameters.

% INPUTS:
% moneyness    -      log-moneyness of the European call options (vector)
% F0           -      forward price at value date (scalar) 
% B            -      discount factor from the maturity to the value date (scalar)
% t            -      maturity date (yearfrac with Act/365) (scalar)
% alpha        -      ATS index of stability
% ATS_params, struct containing:      
%    beta      -      scaling parameter of ATS variance of jumps
%    delta     -      scaling parameter of ATS skew parameter
%    k         -      ATS constant part of the variance of jumps
%    eta       -      ATS constant part of the skew parameter
%    sigma     -      ATS contant diffusion parameter

% OUTPUT:
% prices       -    European call option prices (vector)


%% ATS time-dependent parameters with power-law scaling

k_t   = ATS_params.k   * t^(ATS_params.beta);                              % ATS skew parameter valued in t
eta_t = ATS_params.eta * t^(ATS_params.delta);                             % ATS variance of jumps parameter in t
sigma_t = ATS_params.sigma;                                                % ATS diffusion parameter


%% ATS characteristic function 

phi_t = @(u) exp(  t*(1-alpha)/(k_t*alpha) * ( -1i*u* (1- (1+(k_t*eta_t*sigma_t^2)/(1-alpha))^alpha) + ...
                   + 1 - (1 + k_t*(1i*u*(0.5+eta_t)*sigma_t^2 +0.5*(u*sigma_t).^2 )/(1-alpha)).^alpha)   );


%% Numerical integration

f = @(u) phi_t(-u-1i*0.5)./(u.^2+0.25);                                    % integrand as a function handle

Integral = integral( @(u) exp(-1i*u.*moneyness).*f(u), ...                 % integral in the Lewis formula
                           -inf,inf,'arrayvalued', true) ;

prices = F0*B*( 1 - exp(-moneyness.*0.5)/(2*pi).*real(Integral) );         % European call option prices


end

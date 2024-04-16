function X = LewisFFT_ATS(s,t,alpha,ATS_params,M,N_sim,flag)

% The function simulates one increment of a power-law scaling ATS process between 
% times s and t (with s < t), following the Lewis-FFT simulation algorithm.
% Reference:  Azzone, M. and Baviera, R., A fast Monte Carlo scheme for
%             additive processes and option pricing, 2021. [1]

% INPUTS:
% s            -     first  time for the ATS increment simulation (yearfrac)
% t            -     second time for the ATS increment simulation (yearfrac) 
% alpha        -     ATS index of stability
% ATS_params, struct containing:      
%    beta      -     scaling parameter of ATS variance of jumps
%    delta     -     scaling parameter of ATS skew parameter
%    k         -     ATS constant part of the variance of jumps
%    eta       -     ATS constant part of the skew parameter
%    sigma     -     ATS contant diffusion parameter
% M            -     integer s.t N is the number of grid points
% N_sim        -     number of simulations
% flag         -     1 spline interpolation, 0 linear interpolation

% OUTPUT:
%    X         -      ATS increment between time s and time t


%% Default setting

if nargin < 7
    flag = 1;                                                              % Spline interpolation
end


%% ATS time-dependent parameters with power-law scaling

k_s   = ATS_params.k   * s^(ATS_params.beta);                              % ATS skew parameter valued in s
k_t   = ATS_params.k   * t^(ATS_params.beta);                              % ATS skew parameter valued in t

eta_s = ATS_params.eta * s^(ATS_params.delta);                             % ATS variance of jumps parameter in s
eta_t = ATS_params.eta * t^(ATS_params.delta);                             % ATS variance of jumps parameter in t

sigma_t = ATS_params.sigma;                                                % ATS diffusion parameter


%% Analicity strip bound & Shift
 
p_plus_t = (0.5+eta_t)+sqrt( (0.5+eta_t)^2 ...                             % -(p_plus_t +1) is the lower bound of the strip 
             +2*(1-alpha)/(k_t*sigma_t^2) )-1;                             % of regularity of phi_st [1, page 25]

a = (p_plus_t+1)/2;                                                        % positive real contant (default choice in the equity case)
                                                                           % [1, page 6]

%% Assumption 2: variables which define the characteristic function bound [1, page 25]

if s==0
    
    bound = (1-alpha)^(1-alpha)/(alpha*2^alpha)* ...                       % upper bound of b
            (t*sigma_t^(2*alpha)/(k_t^(1-alpha)));                                    
    log_B = (1-alpha)/(alpha)*(t/(k_t)) + 10;                              % lower bound for B                 
        
else
    
    bound = (1-alpha)^(1-alpha)/(alpha*2^alpha)*   ...                     % upper bound of b
            (  t*sigma_t^(2*alpha)/(k_t^(1-alpha)) ...
              -s*sigma_t^(2*alpha)*(k_s^(alpha-1))  );                      
    log_B = (1-alpha)/(alpha)*(t/(k_t)-s/k_s) + 10;                        % lower bound for B
        
end
   
b = bound - eps;                                                           % 0 < b < bound
omega = 2*alpha - eps;                                                     % 0 < omega < 2*alpha


%% ATS characteristic function 

if s==0 
    phi_st = @(u) exp(  t*(1-alpha)/(k_t*alpha) * ( -1i*u* (1- (1+(k_t*eta_t*sigma_t^2)/(1-alpha))^alpha) + ...
                        + 1 - (1 + k_t*(1i*u*(0.5+eta_t)*sigma_t^2 +0.5*(u*sigma_t).^2 )/(1-alpha)).^alpha)   );
else
    phi_st = @(u) exp(  t*(1-alpha)/(k_t*alpha) * ( -1i*u* (1- (1+(k_t*eta_t*sigma_t^2)/(1-alpha))^alpha) +   ...
                        + 1 - (1 + k_t*(1i*u*(0.5+eta_t)*sigma_t^2 +0.5*(u*sigma_t).^2 )/(1-alpha)).^alpha)   ...
                        - (  s*(1-alpha)/(k_s*alpha) * ( -1i*u* (1- (1+(k_s*eta_s*sigma_t^2)/(1-alpha))^alpha) + ...
                        + 1 - (1 + k_s*(1i*u*(0.5+eta_s)*sigma_t^2 +0.5*(u*sigma_t).^2 )/(1-alpha)).^alpha) )  );
end


%% FFT parameters 

N = 2^M;                                                                   % number of points in the Fourier domain grid

h = (pi*(p_plus_t+1)/(b*N^omega))^(1/(omega+1));                           % step size in the Fourier domain grid
u_1 = 0.5*h;                                                               % first node in the Fourier domain grid
u_grid = u_1 + (0:N-1)*h;                                                  % Fourier domain grid (row vector)

gamma = 2*pi/(h*N);                                                        % step size in the CDF domain grid
x_1 = -(N-1)*gamma/2;                                                      % first node in the CDF domain grid
x_grid = x_1 + (0:N-1)*gamma;                                              % CDF domain grid (row vector)


%% CDF domain grid truncation  

idx = (abs(x_grid-5*sqrt(t-s))==min(abs(x_grid-5*sqrt(t-s))));             % find the index of nearest point to 5*sqrt(t-s) 
x_k = x_grid(idx);                                                         % fix the last node in the truncated CDF domain grid              
xk_grid = (-x_k:gamma:x_k);                                                % truncated CDF domain grid


%% CDF approximation via FFT 

f = @(u) phi_st(u-1i*a)./ (1i*u+a);                                        % integrand as a function handle

f_uj = f( u_grid );                                                        % row vector of f (valued in the nodes u_j)
f_j = f_uj.*exp(-1i*x_1*(0:N-1)*h);                                        % row vector of phases f_j 
FFT_k = fft(f_j);                                                          % FFT computed in f_j nodes (row vector)
f_hat = h*exp(-1i*u_1*x_grid).*FFT_k;                                      % numerical Fourier transform of the function f

Integral = real( f_hat );                                                  % integral in the CDF Lewis formula approximation

Integral = interp1(x_grid,Integral,xk_grid,'spline');                      % interpolate on the truncated grid

P_hat = 1 - exp(-a*xk_grid)/pi .* Integral;                                % Lewis formula for the CDF

% figure('Name','CDF','NumberTitle','off');
% plot(xk_grid,P_hat, 'k-', 'LineWidth',2.5)                                 % plot the CDF                                                            
% title('Cumulative Distribution Function')


%% Assumption 2

% u_grid = [u_grid u_grid(end):h:2000];
% A2 = ( log( abs(phi_st(u_grid-1i*a))) < - b*u_grid.^(omega)+log_B );
% 
% figure('Name','Assumption 2','NumberTitle','off');
% plot( u_grid, A2,'k*', 'LineWidth',1)
% xlabel('u')
% title('Assumption 2 check')


%% Numerical inversion of the CDF via interpolation

U = rand(N_sim,1);                                                         % sample a vector U of N_sim uniform r.v. in [0,1]

[P_hat, idx, ~] = unique(P_hat);                                           % eliminate repetitions to perform the interpolation
xk_grid = xk_grid(idx);                                                    % reduce in a consistent way also the CDF domain grid


switch(flag)
    
    case 0                                                                 % linear interpolation
        X = interp1(P_hat,xk_grid',U,'linear','extrap');                   % interpolate on U
        
    case 1                                                                 % spline interpolation                                                            
        X = interp1(P_hat,xk_grid,U,'spline','extrap');                    % interpolate on U                  
        
    otherwise
        disp('Error: flag must be 0 or 1')
        
end



end

function [S,F]= Underlying_value(discounts,dates,Dates,S0,alpha,ATS_params, M, N_sim,flag)

% The function simulates the dynamics of underlying S and the forward F on a
% discrete-time grid made of settlement date (Dates(1)) reset dates only (Dates(2:end)). 

% INPUTS:
% discounts    -     first  time for the ATS increment simulation (yearfrac)
% dates        -     second time for the ATS increment simulation (yearfrac) 
% Dates        -     vector composed of settlement date and reset dates,
%                    respectively
% S0           -     stock price at time zero (i.e Settlement Date)
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
% S            -     simulated path of the underlying on the discrete-time  
%                    grid described by Dates
% F            -     simulated path of the forward on the discrete-time grid 
%                    described by Dates


%% Default setting

if nargin < 8
    flag = 1;                                                              % Spline interpolation
end


%% Parameters

B_res = GetDiscounts(Dates,dates,discounts);                             % discount factors at settlemene date and reset dates

delta=yearfrac(Dates(1:end-1),Dates(2:end),3);                             % fraction of year between two consecuitve dates of Dates

ttm=yearfrac(Dates(1),Dates,3);                                            % time to maturity corresponding to the reset dates

fwd_zrate=-log(B_res(2:end)./B_res(1:end-1))./delta;                       % forward zero rate


F0=S0/B_res(2);                                                            % initial value of the Forward
F=zeros(N_sim,length(ttm));                                                % initialization of the Forward's simulated path matrix
S=zeros(N_sim,length(ttm));                                                % initialization of the underlying's simulated path matrix
F(:,1)=F0*ones(N_sim,1);                                                   % Forward value at settlement date
S(:,1)=S0*ones(N_sim,1);                                                   % Underlying value at settlement date

for i=1:length(ttm)-2

    rng(i)                                                      

    f_t_s = LewisFFT_ATS(ttm(i),ttm(i+1),alpha,ATS_params, M, N_sim,flag); % ATS increment simulation between two reset dates: t_(i+1) and t_i
    
    S(:,i+1)=F(:,i).*exp(f_t_s);                                           % Underlying value at t_(i+1)                                          
    F(:,i+1)=S(:,i+1).*exp(fwd_zrate(i)*delta(i));                         % Forward value at t_(i+1)

end

index = length(ttm)-1;


rng(index)                                                      

f_t_s = LewisFFT_ATS(ttm(index),ttm(index+1),alpha,ATS_params, M, N_sim,flag); % ATS increment simulation between two reset dates: t_(i+1) and t_i

S(:,index+1)=F(:,index).*exp(f_t_s); 


end

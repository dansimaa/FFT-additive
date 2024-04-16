
% CASE STUDY 01
% Pricing European call option using ATS process increment simulations with
% Lewis-FFT algorithm.


%% Data

SetDate = datenum('19-Feb-2008');                                          % settlement date                                                                                 
Maturity = dateMoveVec(SetDate,'m',1,'MF',eurCalendar);                    % maturity date
B = GetDiscounts(Maturity,dates,discounts);                                % discount factor at maturity date

s = yearfrac(SetDate,SetDate, Act_365);                                    % first  time for the ATS increment simulation
t = yearfrac(SetDate,Maturity,Act_365);                                    % second time for the ATS increment simulation

S0 = 100;                                                                  % stock price on the settlement date
F0 = S0/B;                                                                 % forward price on the settlement date
moneyness = sqrt(t)*linspace(-0.2,0.2,30)';                                % log-moneyness of the 30 Eu call options

%% ATS parameters 

alpha = 3/4;                                                               % ATS index of stability 

ATS_params = struct( 'beta'  ,  1   , ...                                  % scaling parameter of ATS variance of jumps
                     'delta' , -0.5 , ...                                  % scaling parameter of ATS skew parameter
                     'k'     ,  1   , ...                                  % ATS constant part of the variance of jumps
                     'eta'   ,  1   , ...                                  % ATS constant part of the skew parameter
                     'sigma' ,  0.2 );                                     % ATS contant diffusion parameter
                
%% Assumption 1

Assumption1

%% Simulation of the ATS increment

rng(1)                                                                     % set seed
M = 12;                                                                    % integer to define the number of grid points in FFT
N_sim = 10^7;                                                              % number of simulations
flag=1;                                                                    % Lewis-FFT-S algorithm

tic
f = LewisFFT_ATS(s,t,alpha,ATS_params,M,N_sim,flag);                       % ATS increments between s and t
toc
% mean(f)                                                                  % seed(1) -0.001653993163524


%% Pricing

[prices_FFT,CI_FFT] = EuropeanOptionMC_ATS(f, moneyness, F0, B);           % Lewis-FFT-S simulation method prices

prices = LewisFormula_ATS(moneyness, F0, B, t, alpha, ATS_params);         % closed Lewis formula prices


%% Errors

RMSE_error = sqrt(mean((prices_FFT-prices).^2));                           % root mean square error                         

Maximum_error = max(abs(prices_FFT-prices));                               % maximum error                                 

MAPE = 100*mean(abs(prices_FFT-prices)./(prices));                         % mean average percentage error                           


%% Lewis-FFT-S (Spline) vs Lewis-FFT-L (Linear): errors and computational times varying M

M = 7:13;                                                                  % range for the grid parameters

f_st = zeros(N_sim,1);                                                     % initialize the increment vector

prices_FFT_S = zeros(length(moneyness),1);                                 % initialize the price vector (spline) 
prices_FFT_L = zeros(length(moneyness),1);                                 % initialize the price vector (linear)

error_S = zeros(length(M),3);                                              % initialize the error matrix for the spline interpolation
error_L = zeros(length(M),3);                                              % initialize the error matrix for the linear interpolation

Comp_time_S = zeros(length(M),1);                                          % initialize the computational time vector (spline) 
Comp_time_L = zeros(length(M),1);                                          % initialize the computational time vector (linear) 

for ii = 1:length(M)      

    rng(1)
    tic
    f_st = LewisFFT_ATS(s,t,alpha,ATS_params,M(ii),N_sim,0);               % ATS increment with linear interpolation
    Comp_time_L(ii) = toc;                                                 % computational times (linear)
    prices_FFT_L = EuropeanOptionMC_ATS(f_st, moneyness, F0, B);           % Lewis-FFT simulation method prices

    rng(1)
    tic
    f_st = LewisFFT_ATS(s,t,alpha,ATS_params,M(ii),N_sim,1);               % ATS increment with spline interpolation
    Comp_time_S(ii) = toc;                                                 % computational times (spline)
    prices_FFT_S = EuropeanOptionMC_ATS(f_st, moneyness, F0, B);           % Lewis-FFT-S simulation method prices

    error_S(ii,1) = sqrt( mean((prices_FFT_S-prices).^2) );                % root mean square error (spline) 
    error_S(ii,2) = max( abs(prices_FFT_S-prices) );                       % maximum error (spline)
    error_S(ii,3) = 100*mean( abs(prices_FFT_S-prices)./(prices) );        % mean average percentage error (spline)        

    error_L(ii,1) = sqrt(mean((prices_FFT_L-prices).^2));                  % root mean square error (linear)                 
    error_L(ii,2) = max(abs(prices_FFT_L-prices));                         % maximum error (linear)
    error_L(ii,3) = 100*mean(abs(prices_FFT_L-prices)./(prices));          % mean average percentage error (linear) 
    
end


%% Plot the errors in logarithmic scale

for ii = 1:3
    
    figure
    plot(M,log10(error_S(:,ii)),'r*-','LineWidth',2);
    hold on 
    grid on
    plot(M,log10(error_L(:,ii)),'gd-','LineWidth',2);
    legend('spline','linear','FontSize',FntSz)
    xlabel('M','FontSize',FntSz)
    
    if ii==1
        ylabel('log_{10} (RMSE)','FontSize',FntSz)
    end
    if ii==2
        ylabel('log_{10} (MAX error)','FontSize',FntSz)
    end
    if ii==3
        ylabel('log_{10} (MAPE)','FontSize',FntSz)
    end
    
end

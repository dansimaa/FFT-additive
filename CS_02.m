
% CASE STUDY 02
% Fractional Fourier transform (FRFT)


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
                
 
%% Simulation of the ATS increment using FRFT with the same grid as FFT

rng(1)                                                                     % set the seed
M = 12;                                                                    % integer to define the number of grid points in FFT
N_sim = 10^7;                                                              % number of simulations
flag=1;                                                                    % Lewis-FFT-S algorithm

FRFT_alpha = 1/(2^M);                                                      % FRFT parameter to set up the grid

tic
f = LewisFRFT_ATS(s,t,alpha,ATS_params,M,N_sim,FRFT_alpha,flag);           % ATS increments between s and t                   
toc
% mean(f)                                                                  % seed(1) -0.001653993163524
 

%% Pricing 

prices_FRFT = EuropeanOptionMC_ATS(f, moneyness, F0, B);                   % Lewis-FRFT-S simulation method prices

prices = LewisFormula_ATS(moneyness, F0, B, t, alpha, ATS_params);         % closed Lewis formula prices


%% Errors

RMSE_error = sqrt(mean((prices_FRFT-prices).^2));                          % root mean square error                         

Maximum_error = max(abs(prices_FRFT-prices));                              % maximum error                                 

MAPE = 100*mean(abs(prices_FRFT-prices)./(prices));                        % mean average percentage error                           
  

%% Optimization for FRFT_alpha with different metrics
% Error metrics:
%  01) RMSE
%  02) Max error
%  03) MAPE

rng(1)

f_opt = @(frft_alpha) LewisFRFT_ATS(s,t,alpha,ATS_params, M, N_sim,frft_alpha,flag);

distance_01 = @(frft_alpha) sqrt( mean((EuropeanOptionMC_ATS(f_opt(frft_alpha),moneyness,F0,B)-prices).^2) );
distance_02 = @(frft_alpha) max( abs(EuropeanOptionMC_ATS(f_opt(frft_alpha),moneyness,F0,B)-prices) );
distance_03 = @(frft_alpha) 100*mean( abs(EuropeanOptionMC_ATS(f_opt(frft_alpha),moneyness,F0,B)-prices)./(prices) );

FRFT_alpha_01 = fminsearch(distance_01,FRFT_alpha);                        % 2.563476562500000e-04
FRFT_alpha_02 = fminsearch(distance_02,FRFT_alpha);                        % 2.563476562500000e-04
FRFT_alpha_03 = fminsearch(distance_03,FRFT_alpha);                        % 2.563476562500000e-04


%% Errors using the optimized parameters FRFT_alpha for different metrics

rng(1)

f_01 = LewisFRFT_ATS(s,t,alpha,ATS_params, M, N_sim, FRFT_alpha_01,flag);
f_02 = LewisFRFT_ATS(s,t,alpha,ATS_params, M, N_sim, FRFT_alpha_02,flag);
f_03 = LewisFRFT_ATS(s,t,alpha,ATS_params, M, N_sim, FRFT_alpha_03,flag);

prices_FRFT_01 = EuropeanOptionMC_ATS(f_01,moneyness,F0, B);
prices_FRFT_02 = EuropeanOptionMC_ATS(f_02,moneyness,F0, B);
prices_FRFT_03 = EuropeanOptionMC_ATS(f_03,moneyness,F0, B);


RMSE_error_01 = sqrt(mean((prices_FRFT_01-prices).^2));                    % 2.761074772798380e-04
Maximum_error_01 = max(abs(prices_FRFT_01-prices));                        % 5.112134068285812e-04
MAPE_01 = 100*mean(abs(prices_FRFT_01-prices)./(prices));                  % 0.017269121151212

RMSE_error_02 = sqrt(mean((prices_FRFT_02-prices).^2));                    % 2.761074772798380e-04
Maximum_error_02 = max(abs(prices_FRFT_02-prices));                        % 5.112134068285812e-04
MAPE_02 = 100*mean(abs(prices_FRFT_02-prices)./(prices));                  % 0.017269121151212

RMSE_error_03 = sqrt(mean((prices_FRFT_03-prices).^2));                    % 2.761074772798380e-04
Maximum_error_03 =max(abs(prices_FRFT_03-prices));                         % 5.112134068285812e-04
MAPE_03 = 100*mean(abs(prices_FRFT_03-prices)./(prices));                  % 0.017269121151212

%% Optimization procedure for different M values (Spline interpolation)

M = 7:13;                                                                  % range for the grid parameters

f_st_01 = zeros(N_sim,1);                                                  % initialize the increment simulation (RMSE)
f_st_02 = zeros(N_sim,1);                                                  % initialize the increment simulation (MAXERROR)
f_st_03 = zeros(N_sim,1);                                                  % initialize the increment simulation (MAPE)

prices_FRFT_01 = zeros(length(moneyness),length(M));                       % initialize the price matrix (RMSE)
prices_FRFT_02 = zeros(length(moneyness),length(M));                       % initialize the price matrix (MAXERROR)
prices_FRFT_03 = zeros(length(moneyness),length(M));                       % initialize the price matrix (MAPE)

error_matrix_01 = zeros(length(M),3);                                      % initialize the error matrix (RMSE)
error_matrix_02 = zeros(length(M),3);                                      % initialize the error matrix (MAXERROR)
error_matrix_03 = zeros(length(M),3);                                      % initialize the error matrix (MAPE)

Comp_time_01 = zeros(length(M),1);                                         % initialize the computational times vector (RMSE)
Comp_time_02 = zeros(length(M),1);                                         % initialize the computational times vector (MAXERROR)   
Comp_time_03 = zeros(length(M),1);                                         % initialize the computational times vector (MAPE)

FRFT_alpha = zeros(length(M),1);                                           % initialize the alpha vector
FRFT_alpha_01 = zeros(length(M),1);                                        % initialize the optimal alpha vector (RMSE)
FRFT_alpha_02 = zeros(length(M),1);                                        % initialize the optimal alpha vector (MAXERROR)
FRFT_alpha_03 = zeros(length(M),1);                                        % initialize the optimal alpha vector (MAPE)                                      
tic 

for ii = 1:length(M)     

    rng(1)
    FRFT_alpha(ii) = 1/(2^M(ii)); 
   
    f_opt = @(frft_alpha) LewisFRFT_ATS(s,t,alpha,ATS_params, M(ii), N_sim,frft_alpha,flag);
    
    distance_01 = @(frft_alpha) sqrt( mean((EuropeanOptionMC_ATS(f_opt(frft_alpha),moneyness,F0, B)-prices).^2) );
    distance_02 = @(frft_alpha) max( abs(EuropeanOptionMC_ATS(f_opt(frft_alpha),moneyness,F0,B)-prices) );
    distance_03 = @(frft_alpha) 100*mean( abs(EuropeanOptionMC_ATS(f_opt(frft_alpha),moneyness,F0,B)-prices)./(prices) );
    
    FRFT_alpha_01(ii) = fminsearch(distance_01,FRFT_alpha(ii));             
    FRFT_alpha_02(ii) = fminsearch(distance_02,FRFT_alpha(ii));             
    FRFT_alpha_03(ii) = fminsearch(distance_03,FRFT_alpha(ii));             

    rng(1)
    tic
    f_st_01 = LewisFRFT_ATS(s,t,alpha,ATS_params, M(ii), N_sim, FRFT_alpha_01(ii),flag);
    Comp_time_01(ii) = toc;
    prices_FRFT_01(:,ii) = EuropeanOptionMC_ATS(f_st_01, moneyness, F0, B);         

    error_matrix_01(ii,1) = sqrt( mean((prices_FRFT_01(:,ii)-prices).^2) );                
    error_matrix_01(ii,2) = max( abs(prices_FRFT_01(:,ii)-prices) );                       
    error_matrix_01(ii,3) = 100*mean( abs(prices_FRFT_01(:,ii)-prices)./(prices) );            
    
    rng(1)
    tic
    f_st_02 = LewisFRFT_ATS(s,t,alpha,ATS_params, M(ii), N_sim, FRFT_alpha_02(ii),flag);               
    Comp_time_02(ii) = toc;
    prices_FRFT_02(:,ii) = EuropeanOptionMC_ATS(f_st_02, moneyness, F0, B);           

    error_matrix_02(ii,1) = sqrt( mean((prices_FRFT_02(:,ii)-prices).^2) );               
    error_matrix_02(ii,2) = max( abs(prices_FRFT_02(:,ii)-prices) );                      
    error_matrix_02(ii,3) = 100*mean( abs(prices_FRFT_02(:,ii)-prices)./(prices) ); 
    
    rng(1)
    tic
    f_st_03 = LewisFRFT_ATS(s,t,alpha,ATS_params, M(ii), N_sim, FRFT_alpha_03(ii),flag);               
    Comp_time_03(ii) = toc;
    prices_FRFT_03(:,ii) = EuropeanOptionMC_ATS(f_st_03, moneyness, F0, B);           

    error_matrix_03(ii,1) = sqrt( mean((prices_FRFT_03(:,ii)-prices).^2) );                
    error_matrix_03(ii,2) = max( abs(prices_FRFT_03(:,ii)-prices) );                       
    error_matrix_03(ii,3) = 100*mean( abs(prices_FRFT_03(:,ii)-prices)./(prices) );              
    
end

tempo = toc


%% Plot the errors in logarithmic scale

for ii = 1:3
    figure
    plot(M,log10(error_matrix_01(:,ii)),'r*-','LineWidth',2);
    hold on
    grid on
    plot(M,log10(error_matrix_02(:,ii)),'gd-','LineWidth',2);
    plot(M,log10(error_matrix_03(:,ii)),'bo-','LineWidth',2);
    legend('Minimizing RMSE','Minimizing maximum error','Minimizing MAPE','Location','northwest','FontSize',2*FntSz)
    xlabel('M','FontSize',FntSz)
    
    if ii==1
        ylabel('log_{10} (RMSE)','FontSize',2*FntSz)
    end

    if ii==2
        ylabel('log_{10} (MAX error)','FontSize',2*FntSz)
    end

    if ii==3
        ylabel('log_{10} (MAPE)','FontSize',2*FntSz)
    end

end

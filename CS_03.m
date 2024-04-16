
% CASE STUDY 03
% Pricing a Phoenix option 


%% Data
 
SetDate = datenum('19_Feb-2008');                                          % settlement date
Dates = dateMoveVec(SetDate, 'y',(0:6) , 'MF', eurCalendar);               % reset dates
Discount = GetDiscounts(Dates,dates,discounts);                            % discount factors at the reset dates

N_sim = 1e7;                                                               % number of simulations
M = 12;                                                                    % integer to define the number of grid points in FFT
flag = 1;                                                                  % spline interpolation
S0 = 100;                                                                  % stock price on the settlement date


%% ATS parameters

alpha = 3/4;                                                               % ATS index of stability 

ATS_params = struct( 'beta'  ,  1   , ...                                  % scaling parameter of ATS variance of jumps
                     'delta' , -0.5 , ...                                  % scaling parameter of ATS skew parameter
                     'k'     ,  1   , ...                                  % ATS constant part of the variance of jumps
                     'eta'   ,  1   , ...                                  % ATS constant part of the skew parameter
                     'sigma' ,  0.2 );                                     % ATS contant diffusion parameter
                

%% Simulate the underlying value dynamics at the reset dates

tic
[S,F] = Underlying_value(discounts,dates,Dates,S0,alpha,ATS_params, M, N_sim,flag);
Time_ATS = toc;

%% Plot some underlying path realizations 

x=1:7;

figure
hold on
plot(x,S(randi([1 N_sim],3000,1),:),'r-')
grid on 
ylabel('Underlying values','FontSize',2*FntSz)
xlabel('Monitoring dates','FontSize',2*FntSz)


%% Pricing a Phoenix option

k1 = 100;                                                                  % strike price which defines the payoff
k2 = 120;                                                                  % strike price which defines the payoff


[price,price_CI] = PhoenixPricing(S,S0,k1,k2,Discount, N_sim);


%% Black model Phoenix option pricing

DFs = GetDiscounts(Dates(2:end),dates,discounts);                          % discount factors at reset dates
sigma = ATS_params.sigma;                                                  % Black volatility

[price_Black,CI_Black] = PhoenixBlackPricing(SetDate,Dates(2:end),DFs,S0,sigma,k1,k2,N_sim);


%% Computational times comparison

tic
Stock = BlackSimulation(SetDate, Dates(2:end),(S0./DFs), sigma, N_sim);
Time_Black = toc;

disp('ATS underlying path simulation time:')
fprintf(' %.6f   \n', Time_ATS) 
disp('')
disp('Black underlying path simulation time:')
fprintf(' %.6f   \n', Time_Black) 


%% Display

disp('ATS pricing:')
disp([sprintf(' %.6f   ', price_CI(1)) sprintf('%.6f   ', price) sprintf('%.6f', price_CI(2))])
disp('')
disp('Black pricing:')
disp([sprintf(' %.6f   ', CI_Black(1)) sprintf('%.6f   ', price_Black) sprintf('%.6f', CI_Black(2))])

CI_length= abs(price_CI(1)-price_CI(2)) ;


%% Pricing: value and 95% level confidence interval, varying N_sim

% vector of the number of simulations to consider
N_sim_vec = [1e5 2e5 5e5 1e6 2e6 5e6 1e7 2e7];

% Initializing the vectors
TimeVector = zeros(length(N_sim_vec),1);
priceCI_length = zeros(length(N_sim_vec),1);
price_CI = zeros(2,length(N_sim_vec));

% Price and its CI varying N_sim_vec
for i=1:length(N_sim_vec)
    tic
    [S,F]= Underlying_value(discounts,dates,Dates,S0,alpha,ATS_params, M, N_sim_vec(i),1);
    [price,price_CI] = PhoenixPricing(S,S0,k1,k2,Discount, N_sim_vec(i));
    price_CI=price_CI';
    TimeVector(i) = toc;
    priceCI_length(i) = price_CI(2)-price_CI(1);
end


% plot
figure;
yyaxis left
plot(N_sim_vec,priceCI_length,'r*-','LineWidth',2);
ylabel('Length of CI','FontSize',FntSz)
hold on
yyaxis right
plot(N_sim_vec,TimeVector,'bo-', 'LineWidth',2);
ylabel('Elapsed Time','FontSize',FntSz)
xlabel('Number of simulations','FontSize',FntSz)
xticks([1e5 1e6 5e6 1e7 2e7])
legend ('|IC|', 'Time','FontSize',FntSz)
grid on
hold on
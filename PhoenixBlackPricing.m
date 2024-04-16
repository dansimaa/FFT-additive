function [price,CI] = PhoenixBlackPricing(SetDate,Dates,DFs,S0,sigma,K_1,K_2, N)

% The function computes the Phoenix option's price and 95 % IC via a MC
% simulation under the Black model for the equity dynamics. 

% INPUTS:
% setDate       -       settlement dates
% Dates         -       reset dates
% S0            -       equity value at settlement date
% sigma         -       Black volatility
% K_1           -       strike price phoenix
% K_2           -       strike price phoenix
% N             -       number of simulations

% OUTPUTS:
% price         -       simulated Phoenix option price
% price_CI      -       95% IC for the price



%% Forward prices and underlying value prices simulations

Fwd = S0./DFs;                                                             % Forward prices at settlement date
S = BlackSimulation(SetDate, Dates, Fwd, sigma, N);                        % matrix Nxlength(Dates) of simulated equity trajectories


%% Auxiliary matrix for the payoff definition 

I = zeros(size(S));
I(S>K_2) = 1;
I = cumsum(I,2);                                                           % cumulative sum of the rows
I = cumsum(I,2);                                                           % cumulative sum of the rows
I(I>1) = 0;                                                                % only the first element from left is different from zero
index = find(I);                                                           % save the indeces of first elements above K2 for each realization

Zero = ones(N,1)-sum(I,2);                                                 % zero where the sum of the rows in zero
Zero = Zero .* ones(size(S)) ;                                             % zero where the sum of the rows in zero

I = ones(size(S))-cumsum(I,2) - Zero;                                      % fill with one the elements on the left
I(index) = 1;                                                              % fill with one the date which is above K_2


%% Payoff and 95% CI

Payoff = sum( DFs' .* max(S-K_1,0) .* I, 2) / S0 ;

[price, ~, CI] = normfit(Payoff);


end















function [prices,CI] = EuropeanOptionMC_ATS(f, moneyness, F0, B)

% The function computes the prices of European call options, given as input their log-moneyness
% and a simulated increment vector for the log-returns of the forward. The simulation is perfomed
% using the Lewis-FFT algorithm for an ATS process with power-law scaling parameters.

% INPUTS:
% f           -    ATS increment between the value date and the maturity (vector)
% moneyness   -    log-moneyness of the call options (vector)
% F0          -    forward price at value date (scalar)
% B           -    discount factor from the maturity to the value date (scalar)

% OUTPUT:
% prices      -    European call option prices (vector)


%% Pricing 

prices = zeros(length(moneyness),1);                                       % initialize the output vector 

CI = zeros(length(moneyness),2);                                           % initialize the output vector 

for j = 1:length(moneyness)                                                % cicle on the moneyness
    
    w_cc = min( exp(f), exp(-moneyness(j)) );                              % vector of payoffs of a covered call 
    [payoff_cc,~,CI_cc] = normfit(w_cc);                                   % mean of the payoffs of a covered call
    prices(j) = F0*B*(1-payoff_cc);                                        % European call option price for moneyness(ii)
    CI_inv = F0*B*(1-CI_cc');                                             
    CI(j,1) = CI_inv(2);                                                   % CI 95%
    CI(j,2) = CI_inv(1);                                                   % CI 95%
    
end


end
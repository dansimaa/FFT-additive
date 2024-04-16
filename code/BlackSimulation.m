function [S] = BlackSimulation(setDate, ResetDates, Fwd, sigma, N)
rng(1)
% Function which simulates the equity values at input reset dates underl the
% Black 76 model.

% INPUTS:
% setDate       -       settlement dates
% ResetDates    -       reset dates
% Fwd           -       forward prices at settlement date with the reset dates as maturity
% sigma         -       Black volatility
% N             -       number of simulations

% OUTPUTS:
% S             -       equity price simulation at reset dates


%% Define the variables 

Expiries = yearfrac(setDate, ResetDates, 3);                               % Act/365 day count convention

f = zeros(N,1);                                                            % initialize the fwd log-return vector

W = randn(N,1,length(ResetDates));                                         % sample st.n. random vectors of length N

S = zeros(N,length(ResetDates));                                           % initialize the equity price matrix                                  


%% Price simulation 

for ii = 1:length(ResetDates)
    
    if ii==1
        dt = Expiries(ii);                                                          
    else
        dt = (Expiries(ii)-Expiries(ii-1));
    end

    f(:) = f(:) - 0.5 * dt * sigma.^2 + sigma .* W(:,1,ii) * sqrt(dt);     % Black dynamics for the fwd log-returns
    
    S(:,ii) = Fwd(ii) * exp(f(:));                                         % equity value at reset date t_ii
    
end



end

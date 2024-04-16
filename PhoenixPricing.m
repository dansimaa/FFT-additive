function [price,price_CI] = PhoenixPricing(S,S0,k1,k2,Discount, N_sim)
% the function computes the Phoenix option's price and 95 % IC

% INPUTS:

% S : underlying values at settlement and  reset dates
% S0: stock price at t0
% k1: strike (use it in the payoff)
% k2: strike (use it to get when the derivative is called the first time)
% Discount: discount factors at settlement date and reset dates
% N_sim: number of simulations

% OUTPUT:
% price: price (as sample mean) of the option
% price_CI: 95% IC for the price

% For each row, keep that underlying value S just if it is smaller or equal
% to k2 OR it is the first underlying of the row that is bigger than k2

% Preallocating memory for vector index 
index(:,1)=ones(N_sim,1);

% Preallocating memory for vector index 
index_sum=zeros(N_sim,1);

% Preallocating memory for vector index 
payoff=zeros(N_sim,length(Discount)-1);

% Considering the underlying value S for each time interval
for i=1:length(Discount)-1

    % Considering for each row the value of S only if there isn't been a value of S bigger than k2
    index(:,i+1)=index(:,i).*(S(:,i+1)<=k2);

    % Considering for each row the discount factor for the date when S is greater than k2
    index_sum=index_sum+index(:,i);

    % Computing the payoff for each row and for each date (t1,..,t6) only if the value of S is less than k2
    payoff(:,i)=(max((S(:,i+1)-k1),0)).*(index(:,i));

end

% Considering when the underlying is never higher than k2 and so the Phoenix option is never exercised
no_payoff=index(:,end)==1;

% Dividing the payoff matrix by S0 (not computed before since S0 is constant)
payoff=payoff./S0;

% Summing all the payoff for the rows
payoff=sum(payoff,2);

% Considering the payoff only when the Phoenix option is exercised
payoff(no_payoff)=0;

% Discounting the payoff with respect to the first date when S0 is higer than k2 for each row
payoff_attuali=payoff.*Discount(index_sum+1);

% Computing the price and the Confidence interval (95%) for the Phoenix option 
[price,~,price_CI] = normfit(payoff_attuali);

end
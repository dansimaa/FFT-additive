function [dates,discounts]=bootstrap(datesSet,ratesSet)

% function that returns the bootstrap curve made from a dataset

% INPUT   -->   datesSet   -   struct of dates
%                              settlement 
%                              depos
%                              futures
%                              swaps
%               ratesSet   -   struct of rates
%                              depos
%                              futures
%                              swaps

% OUTPUT  -->   dates      -   vector of dates to plot
%               discounts  -   vector of discounts to plot


%% Initialization

dates=datesSet.settlement;                                                                      % I initialize the output          
discounts=1;                                                                                    % I initialize the output  


%% Short of the Discounts Curve (up to 3 months)

Basis=2;                                                                                        % Act/360 convention

datesSet.depos=datesSet.depos(datesSet.depos<=datesSet.futures(1,1));                           % I just take the depos before the first future
YearFrac=yearfrac(dates,datesSet.depos,Basis);                                                  % Year fraction in the Basis convention
ratesSet.depos=ratesSet.depos(datesSet.depos<=datesSet.futures(1,1),:);                         % I just take the depos before the first future

r=mean(ratesSet.depos,2);                                                                       % Observed market deposit rates

discountsDepos=1./(1+YearFrac.*r);                                                              % I compute discount for deposits

dates=[dates; datesSet.depos];                                                                  % I add values to my output
discounts=[discounts; discountsDepos];                                                          % I add values to my output


%% Middle of the Discounts Curve (from 3 months to 2 years)

N=7;                                                                                            % Number of futures considered

YearFrac=yearfrac(datesSet.futures(1:N,1),datesSet.futures(1:N,2),Basis);
r=mean(ratesSet.futures(1:N,:),2);                                                              % Observed market futures rates
discountsFutures=1./(1+YearFrac.*r);                                                            % I compute discount for deposits                                             

Basis=3;                                                                                        % Act/365 convention

SettleDates=datesSet.futures(1:N,1);                                                            % I take dates which goes over the 3 month period
ExpiryDates=datesSet.futures(1:N,2);                                                            % I take dates which goes over the 3 month period

dates=[dates;ExpiryDates];                                                                      % I add values to my output
discounts=[discounts;zeros(N,1)];                                                               % I initialize a part of my output
zRatesStart=zeroRates(dates(1:end-N-1),discounts(1:end-N-1))/100;                               % I compute continuously compounded zero rates not in percentage
zRates=[zRatesStart;zeros(N,1)];

for ii=1:N
    
    zRates(end-N+ii)=zeroRates([dates(1);dates(end-N-1+ii)],...
                               [discounts(1);discounts(end-N-1+ii)])/100;                       % I compute continuously compounded zero rates not in percentage
    zRateExtr=interp1(dates(2:end-N-1+ii),zRates(1:end-N+ii),...
                      SettleDates(ii),'linear','extrap');                                       % Interpolate
    discountFuturesSettle=exp(-zRateExtr*yearfrac(dates(1),SettleDates(ii),Basis));             % Find the discount factor B(t_0,t_ii-1)
  
    discounts(end-N+ii)=discountFuturesSettle*discountsFutures(ii);                             % Find the discount factor B(t_0,t_ii) and update values to my output
   
end


%% End of the Discounts Curve (from 2 years)

N=length(ratesSet.swaps);                                                                       % Number of swap dates

S=mean(ratesSet.swaps,2);                                                                       % Observed market swap rates

SwapDates=datesSet.swaps;                                                                       % I take dates which goes over the 3 month period
dates=[dates;SwapDates];                                                                        % I add values to my output (lunghezza 62)

discounts=[discounts;zeros(N,1)];                                                               % I initialize a part of my output

Basis=3;                                                                                        % Act/365 convention

zRates=zeroRates(dates(1:end-N),discounts(1:end-N))/100;                                        % I compute continuously compounded zero rate not in percentage for the first rate
zRateExtr=interp1(dates(2:end-N),zRates(1:end),SwapDates(1),'linear','extrap');                 % Interpolate
discounts(end-N+1)=exp(-zRateExtr*yearfrac(dates(1),SwapDates(1),Basis));                       % Find the discount factor

Basis=6;                                                                                        % 30/360 European convention
YearFrac=yearfrac(dates([1 end-N+1:end-1]),dates(end-N+1:end),Basis);                           % Compute fraction of years between dates

for ii=2:N                                                                                      % Compute others
    discountSwaps=(1-S(ii)*sum(YearFrac(1:ii-1).*discounts(end-N+1:end-N+ii-1)))/...
                  (1+YearFrac(ii)*S(ii));
    discounts(end-N+ii)=discountSwaps;
end

dates=dates([1:end-N end-N+2:end]);                                                             % I remove the first swap
discounts=discounts([1:end-N end-N+2:end]);


end
















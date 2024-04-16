function DF = GetDiscounts(targetDates,dates,discounts)

% function that returns the discount factors in chosen dates by interpolation on zero
% rates exploiting the DFs curve (dates, discounts) obtained from the bootstrap

% INPUT    -->   targetDates   -   dates (expiries) in which I want to know the discount factor
%                dates         -   dates at which I have the discounts factors, the first element is the settlement date
%                                  dates(1) is the Settlement date
%                discounts     -   discount factors obtained by the bootstrap

% OUTPUTS  -->   DF            -   discount factors at the "target expiries"


%% Conventions

Act_365=3;                                                                  % Act/365 convention


%% Dates

targetDates=datenum(targetDates);                                           % Changing format
dates=datenum(dates);                                                       % Changing format


%% Discounts Factors

zRates        = zeroRates(dates,discounts)/100;                             % Zero rates computation
zRates_interp = interp1(dates,[0;zRates],targetDates,'linear','extrap');    % Interpolating the zero rates in the target expiries
YearFrac      = yearfrac(dates(1),targetDates,Act_365);                     % Compute the yearfrac
DF            = exp(-zRates_interp.*YearFrac);                              % Discount factors computation


end
    


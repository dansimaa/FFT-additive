function zRates = zeroRates(dates, discounts)

% function that returns the continuously compounded zero rates

% INPUT   -->   dates       -   vector
%               discounts   -   vector

% OUTPUT  -->   zRates      -   zero rates in percentage


%% Zero Rates Computation

Basis=3;                                            % Act/365 convention

YearFrac=yearfrac(dates(1),dates(2:end),Basis);     % Year fraction in the Basis convention
zRates=-log(discounts(2:end))./YearFrac*100;        % Zero rates in percentage


end

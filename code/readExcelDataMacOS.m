function [dates, rates] = readExcelDataMacOS( filename)
%  Reads data from excel (MacOS version)
%  It reads bid/ask prices and relevant dates
%  All input rates are in % units
%
% INPUTS:
%  filename: excel file name where data are stored
% 
% OUTPUTS:
%  dates: struct with settlementDate, deposDates, futuresDates, swapDates
%  rates: struct with deposRates, futuresRates, swapRates

%% Dates from Excel

%Define cell-to-array function
cell2arr = @(c) reshape([c{:}], size(c));

%Import all content
allcontent = readcell(filename);

%Settlement date
settlement = allcontent{7,5};
%Date conversion
dates.settlement = datenum(settlement);

%Dates relative to depos
date_depositi = cell2arr(allcontent(10:17, 4));
dates.depos = datenum(date_depositi);

%Dates relative to futures: calc start & end
date_futures_read = cell2arr(allcontent(11:19, 17:18));
numberFutures = size(date_futures_read,1);

dates.futures=ones(numberFutures,2);
dates.futures(:,1) = datenum(date_futures_read(:,1));
dates.futures(:,2) = datenum(date_futures_read(:,2));

%Date relative to swaps: expiry dates
date_swaps = cell2arr(allcontent(38:87, 4));
dates.swaps = datenum(date_swaps);

    
%% Rates from Excel (Bids & Asks)

%Depos
tassi_depositi = cell2arr(allcontent(10:17, 5:6));
rates.depos = tassi_depositi / 100;

%Futures
tassi_futures = cell2arr(allcontent(27:35, 5:6));
%Rates from futures
tassi_futures = 100 - tassi_futures;
rates.futures = tassi_futures / 100;

%Swaps
tassi_swaps = cell2arr(allcontent(38:87, 5:6));
rates.swaps = tassi_swaps / 100;

end % readExcelDataOSX
% Run project
% June 14, 2022
% Chiara Agostini
% Davide Cestaro
% Daniel Sima


clear
close all
clc
rng(1);                                                                                 % Fix a seed to get replicable results
warning('off');                                                                         % No warnings
options = optimset('Display','none');                                                   % Suppress function messages
set(0,'DefaultFigureWindowStyle','docked');                                             % Docking figures (more manageble)
format long                                                                             % More representation precision


%% Settings

formatData = 'dd/mm/yyyy';                                                              % Setting dates format 
FntSz = 20;                                                                             % Setting plots fontsize
FntNm = 'Times';                                                                        % Setting plots fontname


%% Data 

% Read market data
if ispc()                                                                               % Windows version
    [datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap.xls', formatData);
else                                                                                    % MacOS version
    [datesSet, ratesSet] = readExcelDataMacOS('MktData_CurveBootstrap.xls');
end

[dates,discounts]=bootstrap(datesSet,ratesSet);                                         % Bootstrap
                  
Act_360   = 2;                                                                         	% Act/360 convention
Act_365   = 3;                                                                         	% Act/365 convention
EU_30_360 = 6;                                                                       	% 30/360 European convention


%% Case study 1

CS_01


%% Case study 2

CS_02


%% Case study 3

CS_03


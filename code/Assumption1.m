
% Assumption 1 has to be verified for an ATS process:
% p_plus and p_minus are non increasing function in t

tt= linspace(0,300,1000)';                                                 % time discretisation parameter             
 
%% ATS parameters

k_t     =  ATS_params.k   * tt.^(ATS_params.beta);                         % ATS skew parameter valued in t

eta_t   =  ATS_params.eta * tt.^(ATS_params.delta);                        % ATS variance of jumps parameter in t

sigma_t =  ATS_params.sigma;                                               % ATS diffusion parameter


%% Sufficient condition for existence of ATS, first condition [paper 1, page 25]

g1 = (0.5 + eta_t) - sqrt( (0.5 + eta_t).^2 + 2*(1- alpha)./(sigma_t^2.* k_t));

g2 = -(0.5 + eta_t) - sqrt( (0.5 + eta_t).^2 + 2*(1- alpha)./(sigma_t^2.* k_t));

%% p_plus and p_minus [paper 1, page 25]

p_plus =  -g2 -1;
p_minus=  -g1;


%% Graphical representation

figure
plot(tt, p_plus,'r-','LineWidth',1)
hold on
plot(tt, p_minus,'b-','LineWidth',1)
legend ('p+', 'p-','FontSize',2*FntSz)
xlabel('Time','FontSize',2*FntSz)
hold off

%% Numerical check

increment1= (p_plus(2:end)- p_plus(1:end-1));                              % difference between two consecutive values of p_plus  
increment2= (p_minus(2:end)- p_minus(1:end-1));                            % difference between two consecutive values of p_minus  

ver1= isempty(find(increment1>0,1));                                       % look for increasing elements within p_plus
ver2= isempty(find(increment2>0,1));                                       % look for increasing elements within p_minus

if ver1*ver2 ==1                                                           % non increasing behaviour for both p_plus and p_minus
    disp('Assumption 1 holds')
else
    disp('Error: Assumption 1 does NOT hold');
end

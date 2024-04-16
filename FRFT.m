function f = FRFT(h,param_a)

% Fractional Fourier Transform (FRFT) of the input vector h.

% INPUTS:
% h          -     input vector 
% param_a    -     FRFT parameter which controls the grid

% OUTPUT 
% f          -     Fractional Fourier transform of h


N = length(h);                                                             % length of the input vector h

e1 = exp(-pi*1i*param_a*(0:(N-1)).^2)';                                    % auxiliary vector 

e2 = exp(pi*1i*param_a*(N:-1:1).^2)';                                      % auxiliary vector 

y = [h.*e1;  zeros(N,1)];                                                  % vector on which is applied FFT

z = [1./e1;  e2];                                                          % vector on which is applied FFT

Dy = fft(y);                                                               % Fast Fourier Transform of the vector y

Dz = fft(z);                                                               % Fast Fourier Transform of the vector y                                                              

D_k = ifft(Dy.*Dz);                                                        % Inverse Fast Fourier Transform

f = e1.*D_k(1:N);                                                          % FRFT of the input vector h


end
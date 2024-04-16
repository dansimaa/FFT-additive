# FFT Monte Carlo for Additive Processes

This repository contains a MATLAB implementation of the Monte Carlo numerical scheme for additive processes proposed by [Azzone and Baviera (2023)](https://doi.org/10.1007/s10287-023-00463-1). 
The code provides a computational framework to simulate these processes and reproduce the results discussed in the referenced publication.
Additionally, an extension of the algorithm is implemented, exploring the utilization of the fractional fast Fourier transform technique as described by [Chourdakis (2005)](https://www.risk.net/journal-computational-finance/2160574/option-pricing-using-fractional-fft).

## Description

Additive processes are becoming the new frontier in equity derivatives as they are able to accurately replicate market-implied volatility term structures. 
This project led to the implementation of a Monte Carlo numerical scheme leveraging on the fast Fourier transform algorithm which allows to deal with the above mentioned processes. 
This simulation technique is a fast and accurate solution for the path-dependent option valuation problem.

## References

1. Azzone, M., Baviera, R. A fast Monte Carlo scheme for additive processes and option pricing. Comput Manag Sci 20, 31 (2023).[https://doi.org/10.1007/s10287-023-00463-1](https://doi.org/10.1007/s10287-023-00463-1)
2. Chourdakis, K. Option Pricing Using the Fractional FFT. Journal of computational finance 8(2), 1-18 (2005).[link](https://www.risk.net/journal-computational-finance/2160574/option-pricing-using-fractional-fft)

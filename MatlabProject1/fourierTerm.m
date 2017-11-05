function fterm = fourierTerm(ck, k, FO, t)
% function fterm = fourierTerm(ck, k, FO, t)
% fourierTerm calculates and returns a single term of the fourier series
% Inputs:
%   ck = fourier coefficients
%   k = which fourier coefficient
%   FO = fundamental frequency
%   t = vector of times in seconds
% Outputs:
%   fterm = vector of same length t for fourier coefficient

% Do the calculation
fterm = ck*exp(1j*2*k*pi*FO*t);

% Testing using the cosine and sine expansion of eulers
% fterm = ck*(cos(2*pi*k*FO*t) + 1j*sin(2*pi*k*FO*t));

return
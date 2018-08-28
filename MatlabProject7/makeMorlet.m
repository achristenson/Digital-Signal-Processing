function psi = makeMorlet(t,s,f0,ncycles)
% function psi = makeMorlet(t,s,f0,ncycles)
% 
% This function takes in the set of inputs and produces a Morlet wavelet
% for the specified inputs using sigma = ncycles/(2*pi*f0)
% 
% Inputs
%   t = input time vector
%   s = scaling factor
%   f0 = frequency of exponential
%   ncycles = number of cycles for gaussian to decay at f0
% Outputs
%   psi = scaled mother wavelet

% Calculate sigma
sigma = ncycles/(2*pi*f0);

% Calculate the wavelet function
psi = exp(1i*2*pi*f0*t/s).*exp(-1*t.^2/(2*s^2*sigma^2));

return
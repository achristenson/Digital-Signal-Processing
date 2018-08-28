function psi = makeSombrero(t,s)
% function psi = makeSombrero(t,s)
% 
% This function creates a mexican hat wavelet with an optional scaling
% factor
% 
% Inputs:
%   t = time vector
%   s = scaling factor
% Outputs:
%   psi = wavelet
% 

if nargin == 2
    t = t/s;
end

psi = (1-t.^2).*exp(-1*t.^2/2);

return
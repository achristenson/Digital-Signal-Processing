function Dout = DCTcompression(signal,percentRetained)
% function Dout = DCTcompression(signal,percentRetained)
% This function takes in a time domain signal, performs a DCT, retaining
% only a desired percentage of the signal strength.

% Take the n point DCT of the function (making sure it's an even power of
% 2)
n = 2^nextpow2(length(signal));
D = dct(signal,n);

% Find where you want to start zeroing out data
dummy = sort(abs(D));
index = round((1-percentRetained/100)*length(dummy)) + 1;

% Zero out the desired data points
Dout = D;
Dout(abs(Dout)<dummy(index)) = 0;

return
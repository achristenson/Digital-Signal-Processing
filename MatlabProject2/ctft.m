function [X,f] = ctft(x,fs)
%CTFT calculates the continuous-time Fourier transform (CTFT) 
% of a signal x(t) which is sampled at a sampling frequency fs
%
% Usage: [X,f] = ctfts(x,fs)
% 
% The vector X contains the (complex) values of the FT evaluate at the 
% frequencies in in the vector f. In doing so it uses relationships between
% the CTFT and the Discrete Fourier Transform

N = length(x);
X = fftshift(fft(x,N))*(2*pi/N);  % do fft; scale; shift so DC is in the middle
f = linspace(-1,1-1/N,N)*(fs/2);  % return frequencies between -Fs/2 and Fs/2

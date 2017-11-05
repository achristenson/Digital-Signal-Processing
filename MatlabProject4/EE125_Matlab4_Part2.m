%% EE 125 Matlab 4 Part 2: Signal Compression
% Section 1: FFT Compression
% Section 2: Storage savings due to compression
% Section 3: Discrete cosine transform

%% FFT Compression

% Step 1
% Loading the train file
load train

% Listening to the train file
% soundsc(y,Fs)

% Step 2
% Take the N pt FFT using the next highest power of 2 from length(y) for N
n = 2^nextpow2(length(y));
f = Fs*(0:(n/2))/n;
Y = fft(y,n);

% Convert to decibels
P = mag2db(abs(Y));

figure(3)
plot(f,P(1:n/2+1))
xlabel('f (Hz)')
ylabel('|P(f)|, dB')
title('Signal Power vs. Frequency')

% The signal has its most powerful components at 750, 900, 1100, 2100,
% 2700, and 3500 Hz. If you look at the euler decomposition of the Fourier
% transform this means the signal is primarily composed of cosines of those
% frequencies

% Step 3
% Print the FFTcompression function
type('FFTcompression.m')

% Step 4
type('SNRoverall.m')

% Step 5
% Compress using FFT
Y100 = FFTcompression(y,100);
Y50 = FFTcompression(y,50);
Y20 = FFTcompression(y,20);
Y10 = FFTcompression(y,10);

% Reconstruct in time domain
recon100 = ifft(Y100);
recon50 = ifft(Y50);
recon20 = ifft(Y20);
recon10 = ifft(Y10);

% Step 5 a
% Listen to the different sounds
% soundsc(recon100,Fs) % Sounds like the normal signal
% soundsc(recon50,Fs) % Sounds pretty good, but not quite the same
% soundsc(recon20,Fs) % Not terribly different from 50% but not amazing
% soundsc(recon10,Fs) % Not great quality but definitely still sounds like
% a train.
% I did some testing in the command window to try to find where I could
% stop recognizing the signal and I found that to be at about 0.25 percent
% retained, which I found really surprising (I was expecting to not be able
% to remove such a high percentage of the original signal), I could even
% kind of recognize it at 0.1%.

% Step 5 b
% Plot the original magnitude spectrum and the 10% spectrum together
P10 = mag2db(abs(Y10));
P10(P10<-10^6) = 0;
figure(4)
plot(f,P(1:n/2+1),'b',f,P10(1:n/2+1),'r')
xlabel('f (Hz)')
ylabel('|P(f)|, dB')
title('Original spectrum and 10% spectrum compared')
legend('Original','10%')
% The unadulterated signal is much more organic because its composed of a
% much larger variety of sinusoids where the 10% sounds fake because it's
% made up of less harmonics.

% Step 5 c
yp = ifft(Y); % Just doing this to preserve the signal length
snr100 = SNRoverall(yp,recon100);
snr50 = SNRoverall(yp,recon50);
snr20 = SNRoverall(yp,recon20);
snr10 = SNRoverall(yp,recon10);

T = table(snr100,snr50,snr20,snr10)

% Step 6
clear
[y,Fs] = audioread('onscreen.wav');

% Compress using FFT
Y100 = FFTcompression(y,100);
Y50 = FFTcompression(y,50);
Y20 = FFTcompression(y,20);
Y10 = FFTcompression(y,10);

% Reconstruct in time domain
recon100 = ifft(Y100);
recon50 = ifft(Y50);
recon20 = ifft(Y20);
recon10 = ifft(Y10);

% Listening to the outputs
% soundsc(recon100,Fs) % Sounds like the normal signal
% soundsc(recon50,Fs) % Sounds pretty good, guys voice sounds shallower
% soundsc(recon20,Fs) % The guy sounds like he's speaking through an
% intercom
% soundsc(recon10,Fs) % The guy sounds even more fake, like a textbook
% example of a digital recording
% This compression works until just about 2 percent when you start to
% barely be able to understand what is being said

% Plot the original magnitude spectrum and the 10% spectrum together
n = 2^nextpow2(length(y));
f = Fs*(0:(n/2))/n;
Y = fft(y,n);
P = mag2db(abs(Y));
P10 = mag2db(abs(Y10));
P10(P10<-10^6) = 0;
figure(4)
plot(f,P(1:n/2+1),'b',f,P10(1:n/2+1),'r')
xlabel('f (Hz)')
ylabel('|P(f)|, dB')
title('Original spectrum and 10% spectrum compared')
legend('Original','10%')
% Looking at the graph it makes sense why you can still understand what he
% is saying even when retaining only 10 percent of the signal power. This
% is because you can see the frequency waveforms for the important things
% he says like "on screen" and the beep or whatever that is.

% Calculating the SNR values
yp = ifft(Y); % Just doing this to preserve the signal length
snr100 = SNRoverall(yp,recon100);
snr50 = SNRoverall(yp,recon50);
snr20 = SNRoverall(yp,recon20);
snr10 = SNRoverall(yp,recon10);

T = table(snr100,snr50,snr20,snr10)

%% Storage savings due to compression
% It makes sense that the Fourier transform is doubling storage space of
% the signal because the FFT returns the information in both the positive
% and negative frequencies. If I were going to save them, since I know the
% original signal will always be real valued, I would save only half of the
% Fourier transform because the negative frequencies are going to be
% symmetric to the positive frequencies about the y axis. If you look at
% the length adjusted yp, then the amount of data is going to be the same
% whether it is stored in the Fourier or time domain (if you take only one
% half of Fourier data). This would equate to the same amount of memory
% usage storing data in the time or frequency domain assuming you make the
% signal to the next power of 2
whos y
whos yp
whos Y
whos Y100
whos Y50
whos Y20
whos Y10

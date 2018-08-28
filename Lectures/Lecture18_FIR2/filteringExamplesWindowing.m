%%
Fs= 1000;
M=127;
fcutoff = 90;
% now, design hanning and hamming windows

bhann=fir1(M,fcutoff/(Fs/2),hanning(M+1)); % get coefficients for a low-pass filter , omega=0.2 cutoff
[hhann,whann]=freqz(bhann,1);
freqz(bhann,1,1000,Fs); title 'Hanning or Hann'
boldify

bhamm=fir1(M,fcutoff/(Fs/2),hamming(M+1)); % get coefficients for a low-pass filter , omega=0.2 cutoff
hhamm =freqz(bhamm,1);
figure,freqz(bhamm,1,1000,Fs); title 'Hamming'
boldify

% for comparison, add a rectangular
figure,plot(1:M+1,bhann,1:M+1,bhamm,'--'),legend('Hann','Hamming')
title('h(n) plots')
boldify

msgbox('notice that the differences between Hann and Hamming are very subtle in the time domain, but huge in the frequency domain')

%% set up a sinewave signal
nsamp=10000;

f0 = 50;

t = 1/Fs*(1:nsamp);
signal=sin(2*pi*f0*t);

finter = 105;
interference = 2*sin(2*pi*finter*t);

data = signal+interference;

% get frequency spectrum of the signal
[X,f]=ctft(data,Fs);
figure,
subplot(311),plot(t,data) 
title('50 Hz signal + interference')
subplot(312),plot(f,abs(X)), ylabel('|X(\omega)|')
subplot(313),plot(f,20*log10(abs(X)+eps))
ylabel('|X(\omega)|, dB')



%% now see if we remove the interferer
% hanning filter
%b=bhann;
b = bhamm;
a=1;

sig = data; 
% sig = signal;
% filter the signal
y=filter(b,a,sig);
[Y,f]=ctft(y,Fs);

figure,

subplot(311), plot(sig)
xlim([-nsamp/5 nsamp*1.5])
title('input signal')
xlabel('Time, samples')
subplot(312), plot(f,abs(Y))
%xlim([0 Fs/2])
title('filtered signal spectrum')

xlabel('Frequency, Hz (Fs=1kHz')
subplot(313), plot(y)
title('output signal')
xlim([-nsamp/5 nsamp*1.5])
xlabel('Time, samples')
boldify


%% noise example

noise = randn(1,1000);

% get frequency spectrum of the noise
[W,f]=ctft(noise,Fs);


% filter the signal
noiseFilt=filter(b,a,noise);
[Wfilt,f]=ctft(noiseFilt,Fs);

figure,

subplot(221), plot(noise)
title('input noise')
xlabel('Time, samples')
subplot(223), plot(f,20*log10(abs(W)+eps))
%xlim([0 Fs/2])
title('input noise spectrum')

subplot(222), plot(noiseFilt)
title('filtered noise')
xlabel('Time, samples')
subplot(224), plot(f,20*log10(abs(Wfilt)+eps))
%xlim([0 Fs/2])
title('filtered noise spectrum')

xlabel('Time, samples')
boldify

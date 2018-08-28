Fs=1000;
nsamp=17000;
dt=1/Fs;

t = dt*[0:(nsamp-1)];

f1=50;

signal = sin(2*pi*f1*t);
noise = 3*randn(size(signal));

data = signal+noise;

figure,plot(t,data);
hold on, plot(t,signal,'r')

%% try different windows
figure
NFFT=32;
for imult=1:9
    NFFT = NFFT*2;
  
    w = hamming(NFFT);
    windowedData = w.'.*data(1:NFFT);
    
    Y = fft(windowedData,NFFT)/NFFT;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    
    % Plot single-sided amplitude spectrum.
    subplot(3,3,imult)
    plot(f,2*abs(Y(1:NFFT/2+1)))
    titl = sprintf('Single-Sided  Spectrum, %d points',NFFT);
    title(titl)
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
	ylim([0 1])
    pause
end

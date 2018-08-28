Fs=1000;
nsamp=17000;
dt=1/Fs;

t = dt*[0:(nsamp-1)]; 

f1=50;
f2=60;

signal1 = sin(2*pi*f1*t);
signal2 = sin(2*pi*f2*t);

data = signal1 + signal2;


%% try different windows

NFFT=16;
Ncalc=10000;
for imult=1:4
    NFFT = NFFT*2;
  
    w = boxcar(NFFT); disp('boxcar window')
 %   w = hanning(NFFT); disp('hanning window')
    windowedData = w.'.*data(1:NFFT);
    
    Y = fft(windowedData,Ncalc)/NFFT;
    f = Fs/2*linspace(0,1,Ncalc/2+1);
	
	fprintf('Spacing of f vector = %0.2f',median(diff(f)));
    
    % Plot single-sided amplitude spectrum.
    plot(f,2*abs(Y(1:Ncalc/2+1)),'.:')
    titl = sprintf('Single-Sided Spectrum, 10^4 FFT points, %d data points',NFFT);
    title(titl)
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
    pause
end

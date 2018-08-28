Fs=1000;
nsamp=100000;
dt=1/Fs;

whiteNoise = randn(1,nsamp);

% filter the noise to make it colored
[b,a]=butter(2,200/(Fs/2));

data = filter(b,a,whiteNoise);
coloredNoise = data(500:end-500); % eliminate any filter transients

whiteNoise = randn(size(coloredNoise));

if 1
    data = coloredNoise; 
else
    data = whiteNoise;
end

t=dt*(0:length(data)-1);
figure,plot(data)
xlabel('Time, samples')
ylabel('Signal value')
grid
boldify



%% try different windows

NFFT=32;
for imult=1:8
    NFFT = NFFT*2;
  
    %w = hamming(NFFT);
    w = boxcar(NFFT);
    windowedData = w.'.*data(1:NFFT);
    
    Y = fft(windowedData,NFFT);
    f = Fs/2*linspace(0,1,NFFT/2+1);
    
    Sxx = 1./NFFT*abs(Y(1:NFFT/2+1)).^2;
    % Plot single-sided amplitude spectrum.
    plot(f,(Sxx))
    titl = sprintf('Periodogram, %d points',NFFT);
    title(titl)
    xlabel('Frequency (Hz)')
    ylabel('|S_{xx}(f)|^2')
    
    sprintf('mean/variance across spectrum: %2.2f, %2.2f',mean(Sxx),var(Sxx))
    pause
end


%% what if we did averaging? Bartlett method

% fix the window length at M = 1024 points;

NFFT=1024;
f = Fs/2*linspace(0,1,NFFT/2+1);
Sxx_sum = zeros(size(f));

% Bartlett method - no window, no overlap
if 1
    numK=30;
    methodType='Bartlett'; 
    overlap=0;
    w = boxcar(NFFT);
    U=1;
else
    numK=30*2;
    methodType='Welch'; 
    overlap=0.5;
    w = hanning(NFFT);  
    U = sum(w.^2) / NFFT; % so U*NFFT = sum(w^2)
end

barVariance=zeros(1,numK);
figure(102)
for ik=1:numK
    
    istart = 1+(ik-1)*NFFT*(1-overlap) ;
    iend = istart+NFFT-1;
    
    indx = istart:iend;    
 
    windowedData = w.' .* data(indx);
    
    Y = fft(windowedData,NFFT);

    % spectrum for this time
    Sxx_ik = 1./(NFFT*U)*abs(Y(1:NFFT/2+1)).^2;
    
    % running sum
    Sxx_sum = Sxx_sum + Sxx_ik;
    
    % Plot single-sided amplitude spectrum.
    plot(f,Sxx_ik,':',f,Sxx_sum./ik)
    titl = sprintf('%s, NFFT=%d, %d averages',methodType,NFFT,ik);
    title(titl)
    xlabel('Frequency (Hz)')
    ylabel('|S_{xx}(f)|^2')
    legend('single','averaged')
    %boldify
    ylim([0 7])
    grid
    
    barVariance(ik) = var(Sxx_sum./ik);
    sprintf('mean/variance across spectrum: %2.2f, %2.2f',mean(Sxx_sum./ik),var(Sxx_sum./ik))
    
    
    pause(.25)
end

%% what if we did averaging? Welch method

% fix the window length at M = 1024 points;

NFFT=1024;
f = Fs/2*linspace(0,1,NFFT/2+1);
Sxx_sum = zeros(size(f));
figure
% Bartlett method - no window, no overlap
if 0
    numK=30;
    methodType='Bartlett'; 
    overlap=0;
    w = boxcar(NFFT);
    U=1;
else
    numK=30*2;
    methodType='Welch'; 
    overlap=0.5;
    w = hanning(NFFT);  
    U = sum(w.^2) / NFFT; % so U*NFFT = sum(w^2)
end

welchVariance=zeros(1,numK);
figure(103)
for ik=1:numK
    
    istart = 1+(ik-1)*NFFT*(1-overlap) ;
    iend = istart+NFFT-1;
    
    indx = istart:iend;    
 
    windowedData = w.' .* data(indx);
    
    Y = fft(windowedData,NFFT);

    % spectrum for this time
    Sxx_ik = 1./(NFFT*U)*abs(Y(1:NFFT/2+1)).^2;
    
    % running sum
    Sxx_sum = Sxx_sum + Sxx_ik;
    
    % Plot single-sided amplitude spectrum.
    plot(f,Sxx_ik,':',f,Sxx_sum./ik)
    titl = sprintf('%s, NFFT=%d, %d averages',methodType,NFFT,ik);
    title(titl)
    xlabel('Frequency (Hz)')
    ylabel('|S_{xx}(f)|^2')
    legend('single','averaged')
    %boldify
    ylim([0 7])
    grid
    
    welchVariance(ik) = var(Sxx_sum./ik);
    sprintf('mean/variance across spectrum: %2.2f, %2.2f',mean(Sxx_sum./ik),var(Sxx_sum./ik))
    
    
    pause(.25)
end
%%

figure,plot(barVariance,'o'), hold on, plot(welchVariance,'r'), boldify
xlabel('# windows averaged')
ylabel('Variance')
legend('Bartlett','Welch, 0.5 OL')

figure,
plot([1:30]*NFFT,barVariance), hold on, plot([.5:.5:30]*NFFT,welchVariance,'r'), boldify
xlabel('# samples processed')
ylabel('Variance')
xlim([0 30*NFFT])
legend('Bartlett','Welch, 0.5 OL')
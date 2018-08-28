%

fs=8192;
tmax =5; 2;%.75%2;  % 5; for 5 sec, we can here several cyles at 8192
t=linspace(0,tmax,fs*tmax);

omega0 = 2*pi*1500;
beta = 2*pi*3000;

x = sin(omega0*t + .5*beta*t.^2);
frq = (omega0 + beta.*t)/(2*pi);
[Xshort,ws] = ctft(x,fs);

figure,
subplot(211),plot(t,x)
subplot(212),plot(ws,abs(Xshort))
%%
figure,
spectrogram(x, hanning(1024/2),1000/2,[],fs)
caxis([-120 -10])
%caxis([-30 -10])
colorbar

%%

% add on sinusoids
ttone2=0.75;
t=linspace(0,ttone2,fs*ttone2);

omega0 = 2*pi*950;
x2 = sin(omega0*t);

ttone3=0.25;
t=linspace(0,ttone3,fs*ttone3);

omega0 = 2*pi*1000;
omega1 = 2*pi*1600;
x3 = sin(omega0*t);% + sin(omega1*t);

% concatenate them all
xbig = [x x2];% x3];
figure,plot(xbig)
%%
scaleFac = 1%4;  % some interesting numbers: 0.5, 1, 2
figure,
spectrogram(xbig,1024./scaleFac,1000./scaleFac,[],fs,'yaxis')
colorbar
ttl = sprintf('Window length is %2.0f',1024./scaleFac);
title(ttl)

% show "just" top 60 dB
cmx = max(get(gca,'Zlim'));
caxis([cmx-60 cmx]),colorbar

boldify
xlabel('Time, sec')
ylabel('Frequency, Hz')
%%
figure,
N = 64*2;
overlap = round(N*.98);  % 98 percent overlap in windows - nice and smooth
%win = hanning(N); 
win = boxcar(N);
spectrogram(xbig,win',overlap,3*N,fs,'yaxis')
title('boxcar window, 128 points')
cmx = max(get(gca,'Zlim'));
rng=90; 
%rng=15;
caxis([cmx-rng cmx]),colorbar
boldify
xlabel('Time, sec')
ylabel('Frequency, Hz')

win = hann(N);
figure
spectrogram(xbig,win',overlap,3*N,fs,'yaxis')
title('Hann window, 128 points')
cmx = max(get(gca,'Zlim'));
caxis([cmx-rng cmx]),colorbar
boldify
xlabel('Time, sec')
ylabel('Frequency, Hz')
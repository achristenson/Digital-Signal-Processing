
%% Ex1) answer for result worked by hand on board
h=[.5-.5*sin(pi/4) .5 0  -.5 -(.5-.5*sin(pi/4))];

freqz(h,1,1000)
boldify
title('hand-designed filter')

%% make your own bandpass filter
N = 64;
Hd = zeros(1,N);  % put N points around the unit circle
deltaW = 2*pi/N % figure out the spacing
omegas = (0:N-1)/N * 2*pi;
% pick the frequencies from pi/4 up to pi/3
ipass = find(omegas>pi/4 & omegas<=pi/2);
Hd(ipass) = 1; 
% get the negative frequencies
ipass2 = find(omegas>=3*pi/2 & omegas<7*pi/4);
Hd(ipass2) = 1; % make a bandpass!!

figure,plot(omegas,Hd)
title('Desired response')
xlabel('\omega, radians/sec')

h = ifft(Hd);  % get the filter in the time domain
figure,stem(h)

hShift = fftshift(h);
freqz(hShift,1,1000)

%% make your own weird filter
N = 64;
Hd = zeros(1,N);  % put N points around the unit circle
deltaW = 2*pi/N % figure out the spacing
omegas = (0:N-1)/N * 2*pi;
% pick the frequencies from pi/4 up to pi/3
ipass = find(omegas>pi/4 & omegas<=pi/2);
Hd(ipass) = sin(omegas(ipass)).^2; % just for fun - could be anything
% get the negative frequencies
ipass2 = find(omegas>=3*pi/2 & omegas<7*pi/4);
Hd(ipass2) = sin(omegas(ipass2)).^2; % just for fun

figure,plot(omegas,Hd)
title('Desired response')
xlabel('\omega, radians/sec')

h = ifft(Hd);  % get the filter in the time domain
% need to take out tiny imaginary part
h = real(h);
figure,stem(h)

hShift = fftshift(h);
freqz(hShift,1,1000)

%% Ex2) lowpass, fir2
f = [0 0.3 0.5 0.6 0.6 1]; 
m = [1 1 1 1 0 0];
b = fir2(30,f,m,boxcar(31));
[h,w] = freqz(b,1,128);
figure
plot(f,m,':o',w/pi,abs(h))
boldify
legend('Ideal','fir2 Designed')
title('Comparison of Frequency Response Magnitudes - Boxcar')


%% Ex2a) lowpass, fir2 - with Hann window this time
f = [0 0.3 0.5 0.6 0.6 1]; 
m = [1 1 1 1 0 0];
b = fir2(30,f,m,hann(31));
[h,w] = freqz(b,1,128);
figure
plot(f,m,':o',w/pi,abs(h))
boldify
legend('Ideal','fir2 Designed')
title('Comparison of Frequency Response Magnitudes - Boxcar')

%% Ex3) lowpass, fir2, but with 'don't care' region
f = [0 0.3 0.5 0.6 0.7 1]; 
m = [1 1 1 1 0 0];
b = fir2(30,f,m,boxcar(31));
[h,w] = freqz(b,1,128);
figure
plot(f,m,':o',w/pi,abs(h))
boldify
legend('Ideal','fir2 Designed')
title('FIR2 with transition region - Boxcar')

%% Ex4) same as above, but using Hann window
f = [0 0.3 0.5 0.6 0.7 1]; m = [1 1 1 1 0 0];
b = fir2(30,f,m,hann(31));
[h,w] = freqz(b,1,128);
figure
plot(f,m,':o',w/pi,abs(h))
boldify
legend('Ideal','fir2 Designed')
title('FIR2 with transition - Hann')


%% arbirary response, fir2, boxcar
f = [0 0.3 0.5 0.6 0.62 1]; m = [.5 .5 1 .9 0 0];
b = fir2(30,f,m,boxcar(31));
[hbx,w] = freqz(b,1,128);
figure
plot(f,m,':o',w/pi,abs(hbx))
boldify
legend('Ideal','fir2 Designed')
title('Comparison of Frequency Response Magnitudes - Boxcar')


%% arbirary response, fir2, hanning
f = [0 0.3 0.5 0.6 0.62 1]; m = [.5 .5 1 .9 0 0];
b = fir2(30,f,m,hanning(31));
[h,w] = freqz(b,1,128);
figure
plot(f,m,':o',w/pi,abs(h),w/pi,abs(hbx))
boldify
legend('Ideal','fir2, hann','fir2, boxcar')
title('Comparison of Frequency Response Magnitudes - Hann')



%% LEAST SQUARES solver: arbirary response, firls, hanning
f = [0 0.3 0.5 0.6 0.62 1]; m = [.5 .5 1 .9 0 0];
b = firls(30,f,m);
[h,w] = freqz(b,1,128);
figure
plot(f,m,':o',w/pi,abs(h),w/pi,abs(hbx))
boldify
legend('Ideal','firls','fir2, boxcar')
title('Comparison of Frequency Response Magnitudes - least squares')

% check: does FIRLS really give lower error at the specified points?

fSpecBx = interp1(w/pi,abs(hbx),f);
fSpecLs = interp1(w/pi,abs(h),f);

freqSamp_errors = nansum((fSpecBx-m).^2)
LSerrors = nansum((fSpecLs-m).^2)

%% arbirary response, firls vs boxcar: easier problem
f = [0 0.3 0.5 0.6 0.8 1]; m = [.5 .5 1 .9 0 0];
b = firls(30,f,m);
[h,w] = freqz(b,1,128);
bbx = fir2(30,f,m);
[hbx,w] = freqz(bbx,1,128);
figure
plot(f,m,':o',w/pi,abs(h),w/pi,abs(hbx))
boldify
legend('Ideal','firls','fir2, boxcar')
title('Comparison of Frequency Response Magnitudes - least squares')



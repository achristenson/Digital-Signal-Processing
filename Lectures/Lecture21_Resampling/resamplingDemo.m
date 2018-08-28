
%% the 'continuous' signal

%note that it's not really possible to do a real CTFT in matlab - all the
%signals are sampled.  What *is* possible is to use a high sample rate, and
%then convert discrete frequency omega to continuous frequency F so signals
%show up in the right place.

Fs=1000;  % pick a high rate
T = 1/Fs; % the sampling rate
t = (0:10000)*T;  % set up a time vector
x = cos(2*pi*10*t) + 2*cos(2*pi*60*t);


[X,f] = ctft(x,Fs);

figure,
subplot(211),plot(t,x)
xlabel('time, sec')
ylabel('x(t)')
title('x(t) = cos(2*pi*10*t) + 2*cos(2*pi*60*t)')
subplot(212), plot(f,abs(X))
xlim([-100 100])
xlabel('Frequency, Hz')
ylabel('|X(F)|')
boldify


%% sample the signal at 200 Hz

% here ASSUME we apply an anti-alias filter before doing sampling.  Because
% our sample rate is 200 Hz, our anti-alias filter will pass everything
% below 100 Hz
% Note: a real anti-alias filter will roll off slowly, so we'd probably
% pick a lower cutoff

Fs=200;  
T = 1/Fs; % the sampling rate
t = (0:8192-1)*T;  % set up a time vector.  length = power of 2, for reasons we'll discuss later
x1 = cos(2*pi*10*t) + 2*cos(2*pi*60*t);  % sample the signal at these times

X1 = fft(x1);  % take the DTFT
X1 = fftshift(X1);  % matlab gives output over the range [0 , 2pi].  For plotting, shift this to [-pi, pi]

omega = linspace(-pi,pi,length(X1));  % for plotting.  The reason we pick this length will come later in the class...
figure,
subplot(211),stem(t,x1)
xlabel('time, sec')
ylabel('x(t)')
title('x1(n)')

subplot(212), plot(omega,abs(X1))
xlabel('\omega, radians/sec')
ylabel('|X1(\omega)|')
title(' Remember, DTFT repeats periodically outside this range')
boldify



%% sample the signal at 100 Hz, with imaginary perfect anti-alias filter
% 
% % here ASSUME we apply an anti-alias filter before doing sampling.  Because
% % our sample rate is 100 Hz, our anti-alias filter will pass everything
% % below 50 Hz
% % thus, we'll imagine that the 60 Hz signal is completely removed
% 
% Fs=100;  
% T = 1/Fs; % the sampling rate
% t = (0:8192-1)*T;  % set up a time vector.  length = power of 2, for reasons we'll discuss later
% x2 = cos(2*pi*10*t);  % sample the signal at these times
% 
% X2 = fft(x2);  % take the DTFT
% X2 = fftshift(X2);  % matlab gives output over the range [0 , 2pi].  For plotting, shift this to [-pi, pi]
% 
% omega = linspace(-pi,pi,length(X2));  % for plotting.  The reason we pick this length will come later in the class...
% figure,
% subplot(211),stem(t,x2)
% xlabel('time, sec')
% ylabel('x(t)')
% title('x2(n) - after analog antialias filter to remove signals > 50 Hz')
% 
% subplot(212), plot(omega,abs(X2))
% xlabel('\omega, radians/sec')
% ylabel('|X2(\omega)|')
% title(' Remember, DTFT repeats periodically outside this range')
% boldify


%% ATTEMPT 1: downsample x1 by a factor of two without doing any filtering

x3_test = x1(1:2:end);


X3_test = fft(x3_test);  % take the DTFT
X3_test = fftshift(X3_test);  % matlab gives output over the range [0 , 2pi].  For plotting, shift this to [-pi, pi]

omega = linspace(-pi,pi,length(X3_test));  % for plotting.  The reason we pick this length will come later in the class...
figure,
subplot(211),stem(t(1:2:end),x3_test)
xlabel('time, sec')
ylabel('samples')
title('x3(n) test - try just taking every 2nd sample of x1(n)')

subplot(212), plot(omega,abs(X3_test))
xlabel('\omega, radians/sec')
ylabel('|X3test(\omega)|')
title('looks ok - but is the high frequency in the right place?')
boldify

%(2*pi*[10 60])/Fs
%(Fs/(2*pi))*2.5

%% ATTEMPT 2: apply a low-pass filter, then downsample x1 by a factor of


% first, design a filter.   the line below gets us the filter's impulse response
h = fir1(12,0.4)

% apply the filter - for FIR filters, we can just convolve with the impulse
% response
x1f = conv(x1,h,'same');

% the 'same' option is just to make the # samples work nicely for plotting:
%  don't think too hard about it.


% now, we can safely grab every 2nd point
x3 = x1f(1:2:end);


X3 = fft(x3);  % take the DTFT
X3 = fftshift(X3);  % matlab gives output over the range [0 , 2pi].  For plotting, shift this to [-pi, pi]

omega = linspace(-pi,pi,length(X3));  % for plotting.  The reason we pick this length will come later in the class...
figure,
subplot(211),stem(t(1:2:end),x3)
xlabel('time, sec')
ylabel('samples')
title('x3(n)- digital filter, then take every 2nd sample of x1(n)')

subplot(212), plot(omega,abs(X3))
xlabel('\omega, radians/sec')
ylabel('|X3(\omega)|')
title('should match X2 pretty well')
boldify

%% Matlab tricks - 'decimate' vs. 'downsample'vs. 'resample'


Xtest = 1:40;

% downample just takes every Nth sample - no filtering
downsample(Xtest,3)

% decimate applies a filter, then takes every Nth sample
decimate(Xtest,3)

% 'resample' lets you resample at P/Q times the original rate
% so to dowsample by 3, using slightly different filtering than decimate
resample(Xtest,1,3)







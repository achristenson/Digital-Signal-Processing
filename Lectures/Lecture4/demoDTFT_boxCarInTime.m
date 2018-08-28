% boxcar in time

n = 256; % number of points
M = 32;  % number of nonzero points  % CHANGE THIS NUMBER - for example, 16, 32, 64, ....
x = zeros(1,n);

midpoint = n/2 - M/2;
x( (1:M) + midpoint) = 1;  % could use midpoint ~=0 to shift time origin 

figure,subplot(211)
t = [1:n]-n/2;
stem(t, x)
axis tight
title('Rect function in (discrete) time')
xlabel('samples')
ylabel('signal')

subplot(212)
X = fftshift(abs(fft(x)));
omega = linspace(-pi,pi,n);
plot(omega,abs(X));
title('Sinc in (radial) frequency')
xlabel('Frequency \omega')
ylabel('|X(\omega)|')
boldify
grid
axis tight

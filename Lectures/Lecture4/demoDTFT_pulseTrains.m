


%% pulse trains

n = 256;

M = 8;  % CHANGE THIS NUMBER - for example, 16, 32, 64, ....
x = zeros(1,n);
x(1:M:end) = 1;

figure,subplot(211)
stem(x)
axis tight
ttl = sprintf('Pulse train in time, M = %d',M);
title(ttl)
xlabel('samples')
ylabel('signal')


subplot(212)
X = fftshift(fft(x));
omega = linspace(-pi,pi,n);
stem(omega,abs(X));
title('Pulse train in (radial) frequency')
xlabel('Frequency \omega')
ylabel('|X(\omega)|')



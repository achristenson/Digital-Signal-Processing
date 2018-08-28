% boxcar in frequency

n = 256;

omegaC = .1*pi; % CHANGE THIS NUMBER - values > 0 and < pi

% figure out nonzero-components
Nc = round(omegaC * (n/2)/pi);  % we use n points for 2 pi radians
X = zeros(1,n);

% 
X(1:Nc) = 1;
X(end-Nc:end) = 1; % because matlab uses [0, 2pi], not [-pi pi]

figure,subplot(122)
omega = linspace(-pi,pi,n);
plot(omega,fftshift(X));
title('Rect function in (radial) frequency')
xlabel('Frequency \omega')
ylabel('|X(\omega)|')
boldify
grid
axis tight

subplot(121),
time = [1:n]-n/2;
stem(time,fftshift(abs(ifft(X))))
axis tight
title('Sinc in (discrete) time')
xlabel('samples')
ylabel('signal')


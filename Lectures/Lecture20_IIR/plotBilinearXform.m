

% first, plot H(s) from Ex 10.3.1
F = linspace(0,4,100);

s = j*2*pi*F;
Ha = (s+0.1)./( (s+.1).^2 + 16);

figure,
subplot(211),plot(2*pi*F,abs(Ha))
xlabel('2*pi*Frequency, Hz')
ylabel('|Ha(F)|')
subplot(212),plot(2*pi*F,unwrap(angle(Ha)))
xlabel('2*pi*Frequency, Hz')
ylabel('phase Ha(F)')


% now look at digital equivalent using bilinear


omega = linspace(0,pi,100); %2*pi/Fs * F;  % digital frequency
z = exp(j*omega);

sApprox = 4*(1-1./z)./(1+1./z);
Hd = (sApprox+0.1)./( (sApprox+.1).^2 + 16);


figure,
subplot(211),plot(omega,abs(Hd))
xlabel('Frequency, rad/sec')
ylabel('|Hz(\omega)|')
subplot(212),plot(omega,unwrap(angle(Hd)))
xlabel('Frequency, rad/sec')
ylabel('phase Hz')
title('Digital IIR design by Bilinear')


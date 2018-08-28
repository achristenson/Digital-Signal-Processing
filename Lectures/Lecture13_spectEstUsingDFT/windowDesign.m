
%  rectangular
bx32=boxcar(32);
[W32,omega32]=freqz(bx32./sum(bx32),1);  % get frequncy response of (normalized) window alone:  W(w)


%  rectangular
bx62=boxcar(62);
[W62,omega62]=freqz(bx62./sum(bx62),1);  % get frequncy response of (normalized) window alone:  W(w)

% triangular window vs. rectangular
bart62=bartlett(62);
Wbart62=freqz(bart62./sum(bart62),1);  % get frequncy response of (normalized) window alone:  W(w)


figure,
plot(bartlett(62))
grid, boldify
title('Bartlett window')

figure,
plot(omega32/pi,20*log10(abs(W32)),omega62/pi,20*log10(abs(W62)), omega62/pi,20*log10(abs(Wbart62)),'r--')
xlabel('\omega/ \pi')
ylabel('|W(\omega)|, dB')
grid
legend('M=32','M=62','Bartlett, 62')
ylim([-150 25])
boldify


%% Hanning vs. rectangular
han62=hanning(62);
Whan62=freqz(han62./sum(han62),1);  % get frequncy response of (normalized) window alone:  W(w)

figure,
plot(han62)
grid, boldify
title('Hanning window')

figure,
plot(omega62/pi,20*log10(abs(Whan62)),omega62/pi,20*log10(abs(W62)),'r--')
xlabel('\omega/ \pi')
ylabel('|W(\omega)|, dB')
grid
legend('Hanning, 62','boxcar, 62')
ylim([-150 25])
boldify

%% Hamming vs. rectangular
ham62=hamming(62);
Wham62=freqz(ham62./sum(ham62),1);  % get frequncy response of (normalized) window alone:  W(w)

figure
plot(ham62)
hold on,
plot(han62,'--')
grid, boldify
legend('Hamming window','Hanning')

figure,
plot(omega62/pi,20*log10(abs(Whan62)),...
omega62/pi,20*log10(abs(Wham62)),'m-.',...
omega62/pi,20*log10(abs(W62)),'r--')
xlabel('\omega/ \pi')
ylabel('|W(\omega)|, dB')
grid
legend('hanning, 62','hamming (M), 62','boxcar')
ylim([-150 25])
boldify


%% Kaiser examples

n=62;

w1=kaiser(n,.5);
w2=kaiser(n,2.5);
w3=kaiser(n,7.5);
npt=length(w1);
figure,
plot(1:npt,w1,1:npt,w2,1:npt,w3)
grid, boldify
title('Kaiser window')
legend('\beta=0.5','\beta=2.5','\beta=7.5')

[Hw1,omega1]=freqz(w1./sum(w1),1);
Hw2 = freqz(w2./sum(w2),1);
Hw3 = freqz(w3./sum(w3),1);


figure,
%subplot(211)
plot(omega1/pi,20*log10(abs(Hw1)),omega1/pi,20*log10(abs(Hw2)),omega1/pi,20*log10(abs(Hw3)))
xlabel('\omega/ \pi')
ylabel('|W(\omega)|, dB')
legend('\beta=0.5','\beta=2.5','\beta=7.5')
%legend('\beta=0.5','\beta=0.2','\beta=0.8')
grid, boldify
title('Kaiser-window response')


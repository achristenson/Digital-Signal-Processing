

% first, plot H(s) from Ex 10.3.1
F = linspace(0,4,100);

s = j*2*pi*F;
Ha = 1./( (s+.1).^2 + 0.9); 

% now look at digital equivalent using bilinear

% interesting values:  10, 2, 1 Hz
Fs_vec = [100 10 5 2];
for k = 1:length(Fs_vec)
    Fs = Fs_vec(k);  % sampling rate, Hz
    T = 1/Fs;
    omega = 2*pi/Fs * F;  % digital frequency
    z = exp(j*omega);
    sApprox = 2/T*(1-1./z)./(1+1./z);
    Hz = 1./( ( sApprox + .1).^2 + 0.9);
    
    figure,
    subplot(211),plot(F,abs(Hz),F,abs(Ha),'.')
    xlabel('Frequency, Hz')
    ylabel('|Hz(\omega)|')
    title('Digital IIR design by Bilinear transform')
    legend('digital','analog')
    
    subplot(212),plot(F,unwrap(angle(Hz)),F,unwrap(angle(Ha)),'.')
    xlabel('Frequency, Hz')
    ylabel('phase Hz')
    title(['sampling rate ' num2str(Fs) ' Hz'])
    
    pause 
end


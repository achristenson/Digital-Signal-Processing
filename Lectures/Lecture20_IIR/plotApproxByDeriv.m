

% first, plot H(s) from Ex 10.3.1
F = linspace(0,1,100);

s = j*2*pi*F;
Ha = 1./( (s+.1).^2 + 0.9); 


% now look at digital equivalent

% interesting values:  10, 2, 1 Hz
Fs_vec = [100 10 5 2 1];
for k = 1:length(Fs_vec)
    Fs = Fs_vec(k);  % sampling rate, Hz
    T = 1/Fs;
    omega = 2*pi/Fs * F;  % digital frequency
    z = exp(j*omega);
    Hz = 1./( ( (1-1./z)/T  + .1).^2 + 0.9);
    
    figure,
    subplot(211),plot(F,abs(Hz),F,abs(Ha),'.')
    xlabel('Frequency, Hz')
    ylabel('|Hz(\omega)|')
    title('Digital IIR design by Approx of Derivatives')
    legend('digital','analog')
    subplot(212),plot(F,unwrap(angle(Hz)),F,unwrap(angle(Ha)),'.')
    xlabel('Frequency, Hz')
    ylabel('phase Hz')
    title(['sampling rate ' num2str(Fs) ' Hz'])
    pause 
end


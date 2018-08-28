

L = 40;
hn = hann(L);

bk = blackman(L);

beta = [.1 .5 1 1.5 2 2.5 3 3.5 4 4.5 5 6  7 9 12 15 ];

figure
ii = 1:L;
for ik=1:length(beta);
    w = kaiser(L,beta(ik));
     

    subplot(311)
    plot(ii,w,ii,hn,ii,bk)
    legend('Kaiser','Hann','Blackman')
    ttl = sprintf('Kaiser window, beta = %2.2f',beta(ik));
    title(ttl)
    boldify
    grid on
    
    subplot(312)
    wnorm = w./sum(w);
    [W,omega] = freqz(wnorm,1,1024);
    plot(omega/pi,20*log10(abs(W)))
    grid on 
    boldify
    ylim([-100 10])
    xlabel('normalized frequency (omega / pi)')
    ylabel('|W(\omega)|')
    
    % now, use this window to design a bandstop filter
    Fs = 1000;
    %bandpassEdges = [245 250];
    bandpassEdges = [200 300];
    b=fir1(L-1,bandpassEdges./(Fs/2),w);
    [H,f] = freqz(b,1,128,Fs);
    subplot(313)
    plot(f,20*log10(abs(H)))
    ylim([-150 25])    
    grid on 
    boldify    
    xlabel('Frequency, Hz')
    ylabel('|H(\omega)|')
    title('Kaiser window used in bandpass')
    
    pause
end
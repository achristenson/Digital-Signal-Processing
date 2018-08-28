M = 101;                           
normF = [0 0.3 0.4 0.6 0.8 1.0]; % transition bands different
%normF = [0 0.3 0.4 0.6 0.7 1.0];  % transition bands the same
amp = [0 0 1 1 0 0];              % desired amplitude in each band

[b2,err2] = firpm(M-1,normF,amp); % optimal filter of length M


figure(1);
[H2,freq]=freqz(b2,[1],512);
subplot(2,1,1); plot(freq*10000/pi,20*log10(abs(H2))); grid;
title('Parks-McClellen, with problem')
xlabel('Hz'); ylabel('dB'); axis([0 10000 -110 30]);
subplot(2,1,2); plot(freq*10000/pi,20*log10(abs(H2))); grid;
xlabel('Hz'); ylabel('dB'); axis([3500 7000 -0.02 0.02]);
boldify

%%
normF = [0 0.3 0.4 0.6 0.7 1.0];  % transition bands the same  % old was .... 0.8 1.0]
amp = [0 0 1 1 0 0];              % desired amplitude in each band

[b2,err2] = firpm(M-1,normF,amp); % optimal filter of length M


figure(2);
[H2,freq]=freqz(b2,[1],512);
subplot(2,1,1); plot(freq*10000/pi,20*log10(abs(H2))); grid;
title('Parks-McClellen, fixed')
xlabel('Hz'); ylabel('dB'); axis([0 10000 -110 30]);
subplot(2,1,2); plot(freq*10000/pi,20*log10(abs(H2))); grid;
xlabel('Hz'); ylabel('dB'); axis([3500 7000 -0.02 0.02]);
boldify


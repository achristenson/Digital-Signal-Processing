function snr = SNRoverall(signal,reconSignal)
% function snr = SNRoverall(signal,reconSignal)
% This function computes the overall SNR of a signal

snr = 10*log(sum(signal.^2./(signal-reconSignal).^2));

return
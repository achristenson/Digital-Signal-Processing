function [ampl,latency] = analyzeSNAP(x,Fs)
% function [ampl,latency] = analyzeSNAP(x,Fs)
% This is a function that takes in the SNAP (Sensory nerve action pulse)
% data (uV) and the sample rate (kHz) for the data and returns the latency
% (ms) and the amplitude (uV) for the set

[MAX,iMax] = max(x);
[MIN,iMin] = min(x);
ampl = (MAX-MIN)/2;
latency = (iMax+iMin)/(2*Fs);

return


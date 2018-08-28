function [data,fs] = demoData

fs=200;
tmax = 1;
t=linspace(0,tmax,fs*tmax);

% set up two tones at 10 and 15 Hz
x1 = sin(2*pi*10*t) + sin(2*pi*20*t);

% set up some bursts
burst = zeros(size(x1));
ion = 50:60;
burst(ion) = sin(2*pi*40*t(ion));

burst2 = zeros(size(x1));
ion = 120:130;
burst2(ion) = sin(2*pi*60*t(ion));

data = x1 + burst + burst2;

return
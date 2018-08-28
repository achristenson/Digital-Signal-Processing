Fs = 3000;

t = linspace(0,2,2*Fs); % set up time from 0-2 sec at desired Fs

y = sin(2*pi*600*t + pi/10);  % 600 Hz sinewave

% then do these:
% sound(y,Fs)
% sound(y(1:2:end),Fs/2)
% sound(y,Fs)
% sound(y(1:3:end),Fs/3)
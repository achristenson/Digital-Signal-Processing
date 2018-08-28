clear all
close all
clc

L = 32;
n = 0:L-1;
x = cos(pi/16*n);  % also compare to pi/15

X = fft(x);
f = linspace(0,2,length(X)+1);
f = f(1:end-1); % above, I used length + 1 to get a redundant point at 2 pi
% in line above, discard the redundant point

figure,stem(f,abs(X));
%ylim([0 11])
ylabel('|X|')
xlabel('\omega/\pi')
%boldify
grid

%%

Xpad=fft(x,L*100);
f = linspace(0,2,length(Xpad)+1);
f = f(1:end-1); % above, I used length + 1 to get a redundant point at 2 pi
% in line above, discard the redundant point
figure,stem(f,abs(Xpad));
ylabel('|X|');
xlabel('\omega/\pi')
%boldify
grid
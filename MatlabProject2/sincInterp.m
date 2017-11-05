function y = sincInterp(x,T,t)
% This function takes in a vector x, period T, and time vector t and
% creates a sinc interpolation y from the data

L = length(t);

y = zeros(1,L);

% Assuming causal signal so can start at n=0
for n = 0:length(x)-1
    d = x(n+1)*sinc(pi*(t-n*T)/T);
    y = y+d;
end

return
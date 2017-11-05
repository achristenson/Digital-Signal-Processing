function vecSum = sumFourierSeries(K,Tp,tau)
% function vecSum = sumFourierSeries(K,Tp,tau)

vecSum = 0;

% set up a time vector t
t = -4:(1/100):4;

% now, loop over the coefficients k, from -K to K
% Inside the loop, do:
%    calculate the value of ck, from pulseTrainDFS.m
%    calculate the fourier term from 'fourierTerm.m'
%    add up the new term to a running sum 'vecSum'
for k =-K:K
    ck = pulseTrainDFS(k,Tp,tau);
    test = fourierTerm(ck,k,Tp,t);
    vecSum = vecSum + test;
end

return
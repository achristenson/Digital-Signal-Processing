function ck = pulseTrainDFS(k,Tp,tau)
% function ck = pulseTrainDFS(Tp,tau)
% implements Fourier series as calculated in example 4.1.1 of P&M, 4th
% edition  
% Inputs:  
%    k, the coefficient index. This function will work if you enter 
%    EITHER a single value for k (i.e., '4') or a vector (i.e.,'-5:5')
%    Tp, the fundamental period of the pulse train (see Fig 4.1.3)
%    tau, the width of the each pulse
%
% ASSUMES the amplitude A = 1

% formula assumes both Tp and tau aren't zero; check
assert(Tp*tau~=0,'Tp and tau must both be nonzero')

F0 = 1/Tp;
% above:if the pulse train is repeating at times Tp (fig 4.1.3), 
% the fundamental period is Tp, so the fundamental freq is 1/Tp

% do the calculation in Eq 4.1.18.  
% Note, this gives the WRONG ANSWER for k=0
ck = tau/Tp *sin(pi*k*F0*tau)./(pi*k*F0*tau);

% go back and correct the k=0 term
ck(k==0) =  tau/Tp;  % eq 4.1.17 in book


return
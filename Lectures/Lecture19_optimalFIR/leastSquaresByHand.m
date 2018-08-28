

% constraint 1: we want a linear phase, causal filter.
% Solution: we'll cheat and design a real, symmetric filter and then will time-shift
% it after our design is finished.  This means that we will enforce
% w(n) = h(n), n>=0 AND
%      = h(-n), n<0
% thus for an odd filter, length M, we will have Mh = (M-1)/2 symmetric
% points on either side of zero
% because of this, W(omega) = h(0) + 2 sum_{n=1}^{n=Mh) h(n) cos(omega n)

% set the filter length
M = 91;
Mh = (M-1)/2; 

% set up a vector of frequencies:
nFreq = 100;
f = linspace(0,1,nFreq);
omega = f*pi;  % zero to pi

% set up the desired response at those frequencies
Hd = zeros(nFreq,1);
Hd(f<0.4) =1;



% set up the matrix that multiplies the h(n) values
A = zeros(nFreq,Mh+1);
A(:,1) = 1;  % multiplies h(0) at each frequency
for ifreq=1:nFreq
    for n = 1:Mh
        A(ifreq,n+1) = 2*cos(omega(ifreq)*n);
    end
end

%% try just solving directly - normal least-squares
hpos = A\Hd;

% For homework, don't do above... instead, use the pseudo inverse
% explicitly!


% now, stick back the other side
hneg = flipud(hpos(2:end));
b = [hneg; hpos];

figure
freqz(b,1)
title('direct solution: = least squares if A is not square')

%% weighted least squares
%weigthing vector
w = ones(size(Hd));
w(f>0.5 & f<0.6)=100;

hposWLS = lscov(A,Hd,w);  % use the built in matlab tool

% now, stick back the other side
hnegWLS = flipud(hposWLS(2:end));
bWLS = [hnegWLS; hposWLS];

bWLS = bWLS./sum(bWLS);
freqz(bWLS,1)
title('*weighted* least squares if A is not square')
%% do the same thing by hand
bigW = diag(sqrt(w));
wA = bigW*A;
wHd = bigW*Hd;

w_hpos = wA\wHd;



% now, stick back the other side
w_hneg = flipud(w_hpos(2:end));
wb = [w_hneg; w_hpos];

figure
freqz(wb,1)
title('direct solution: = weighted least squares if A is not square')








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
f = linspace(0,1,100);
fcut = 0.4;
% cut out frequncies around this
f(f>0.375 & f<0.425) = 0;

nFreq = length(f);

omega = f*pi;  % zero to pi

% set up the desired response at those frequencies


d = zeros(nFreq,1);
d(f<=fcut) =1;



% set up the matrix that multiplies the h(n) values
A = zeros(nFreq,Mh+1);
A(:,1) = 1;  % multiplies h(0) at each frequency
for ifreq=1:nFreq
    for n = 1:Mh
        A(ifreq,n+1) = 2*cos(omega(ifreq)*n);
    end
end



%% try just solving directly - normal least-squares
hpos = A\d;

% prove that this is the pseudoinverse


% now, stick back the other side
hneg = flipud(hpos(2:end));
b = [hneg; hpos];

figure
[hLS,omega] = freqz(b,1);
figure,plot(omega./pi,mag2db(abs(hLS)))

title('direct solution: = least squares if A is not square')

%% weighted least squares
%weigthing vector
 
% single iteration to reduce larger errors
% 1) compute error from above
e = abs(A*hpos - d).^2;


% plot the error

% now, give more weight to regions of high error
w = e;
% but add in a don't care region
%w(f>0.38 & f<0.42) = 0;

figure
plot(f,e,f,w)
xlabel('Frequency')
ylabel('squared error (or weight)')
legend('error','weight')

hposWLS = lscov(A,d,w);

% now, stick back the other side
hnegWLS = flipud(hposWLS(2:end));
bWLS = [hnegWLS; hposWLS];


hReweight = freqz(bWLS,1);

figure,
plot(omega./pi,mag2db(abs(hLS)),omega./pi,mag2db(abs(hReweight)))
legend('LS','reWLS1')
title('*weighted* least squares if A is not square')


%% now, try IRLS
p = 70;
K  = 1.4;
numIter = 20;
birls = IRLS1(A,d,p,K,numIter);

bRWLS = [flipud(birls(2:end)); birls];
hIRLS = freqz(bRWLS);

figure,
plot(omega./pi,mag2db(abs(hLS)),omega./pi,mag2db(abs(hIRLS)))
legend('LS','IRLS,')
title(sprintf('IRLS, %d-norm',p))



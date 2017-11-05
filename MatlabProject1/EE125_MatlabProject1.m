% Problem 1
% Section 1.2
% Initialize the time vector
% Sampling frequency in Hz
fs = 100;
% Creating the time vector from 0 to 1
t = 0:(1/fs):1;

% Plot all of the graphs
figure(1);
% Initialize the variables for the fourier calculations
ck = 1;
FO = 1/3;
% Using a for look for the values of k since they are consecutive and nice
% Note that there will be 8 total graphs
for k=1:3
    test = fourierTerm(ck, k, FO, t);
    subplot(4,2,(2*k-1));
    plot(t,test);
    title(sprintf('Plot(t,test) with k = %d',k));
    subplot(4,2,2*k);
    plot(test);
    title(sprintf('Plot(test) with k = %d',k))
end
% Comparing the graphs for k = FO = 1
k = 1;
FO = 1;
test = fourierTerm(ck,k,FO,t);
subplot(4,2,7);
plot(t,test);
title('Plot(t,test) with k = Fo = 1');
subplot(4,2,8);
plot(test);
title('Plot(test) with k = Fo = 1');

% Discussion
% Modifying the k term in the exponent modifies the frequency of the
% sinusoids that are being plotted. When you plot the fourier term with
% respect to time you ignore the imaginary components of them, where just
% using the plot function includes the imaginary part.

% Section 1.3, part 1
% Create the time vector from 0 to 2 sec, sampled at 100Hz
t = 0:(1/fs):2;

% Set up the coefficients for the fourier term
ck = 1;
k = 1;
FO = 2;

% Call the fourier term function
test = fourierTerm(ck,k,FO,t);

% Plot the resulting data
figure(2);
plotComplexData(t,test);

% Discussion part 1
% The real and imaginary parts represent the real cosine part and the sine
% imaginary part of the Fourier term. The absolute value should be constant
% because there is no time dependent attenuation on the Fourier term and
% the unwrapped phase is a straight line because there is a constant time
% term affecting the phase between the real and the imaginary components of
% the Fourier term. If the phase were wrapped it would be the unit circle.

% Section 1.3, part 2
% Update the coefficients for the fourier term
ck = exp(1j*pi/4);

% Call the fourier term function
test = fourierTerm(ck,k,FO,t);

% Plot the resulting data
figure(3);
plotComplexData(t,test);

% Discussion part 2
% The situation where you have an imaginary Fourier coefficient ck, there
% will be a phase shift introduced to the term. This is evident in the
% plots of the real and imaginary components with respect to time. A
% vertical shift in the phase plot corresponds to a horizontal shift in the
% time domain.

% Section 1.3, part 3
% Update the coefficients for the fourier term
ck = 1;

% Call the fourier term function
test = fourierTerm(ck,k,FO,t) + fourierTerm(ck,k,-1*FO,t);

% Plot the resulting data
figure(4);
plotComplexData(t,test);

% Discussion part 3
% The phases for these Fourier terms are designed to cancel out the 
% imaginary components of each other but leave the real components intact.
% This translates into the euler expansion because sine is an odd function
% and cosine is an even function, so when you make the fundamental
% frequencies opposite of each other and then add the components, the
% cosine parts add and the sine parts become opposites of each other and
% cancel out.

% Section 1.3, part 4
% Update the coefficients for the fourier term
ck = 1*j;

% Call the fourier term function
test = fourierTerm(ck,k,FO,t) + fourierTerm(-1*ck,k,FO,t);

% Plot the resulting data
figure(5);
plotComplexData(t,test);

% Discussion part 4
% Similar to the last situation, in the euler expansion the sine components
% cancel each other out, but additionally the cosine components cancel each
% other, so the magnitude of the signal becomes zero. This is because the
% fourier coefficients make the Fourier terms complex coefficients of each
% other so they are exactly opposite and cancel out perfectly.

% Section 1.3, part 5
% Update the Fourier term and plot
test = test + fourierTerm(2,0,FO,t);
figure(6);
plotComplexData(t,test);

% Discussion part 5
% Adding another fourier term with co = 2 for k = 0 does nothing to the
% fourier term besides add a vertical shift to the real component to it.
% This is basically adding a 0 frequency component to the Fourier term. In
% this case the fundamental frequency is irrelevant because the k=0 term
% forces it to be zero regardless of what FO is.

type('fourierTerm.m');
type('plotComplexData.m');

% Problem 2 Synthesizing other waveforms
% Choosing values for Tp and tau
Tp = 2;
tau = 0.5;

% Compute x(t) for K = 1,2,3,10,50,100,500
K = 1;
X = sumFourierSeries(K,Tp,tau);

figure(7);
plot(real(X));
xlabel('time (centiseconds)');
ylabel('x(t)');
title(sprintf('x(t) vs. time, K = %d',K));

K = 2;
X = sumFourierSeries(K,Tp,tau);

figure(8);
plot(real(X));
xlabel('time (centiseconds)');
ylabel('x(t)');
title(sprintf('x(t) vs. time, K = %d',K));

K = 3;
X = sumFourierSeries(K,Tp,tau);

figure(9);
plot(real(X));
xlabel('time (centiseconds)');
ylabel('x(t)');
title(sprintf('x(t) vs. time, K = %d',K));

K = 10;
X = sumFourierSeries(K,Tp,tau);

figure(10);
plot(real(X));
xlabel('time (centiseconds)');
ylabel('x(t)');
title(sprintf('x(t) vs. time, K = %d',K));

K = 50;
X = sumFourierSeries(K,Tp,tau);

figure(11);
plot(real(X));
xlabel('time (centiseconds)');
ylabel('x(t)');
title(sprintf('x(t) vs. time, K = %d',K));

K = 100;
X = sumFourierSeries(K,Tp,tau);

figure(12);
plot(real(X));
xlabel('time (centiseconds)');
ylabel('x(t)');
title(sprintf('x(t) vs. time, K = %d',K));

K = 500;
X = sumFourierSeries(K,Tp,tau);

figure(13);
plot(real(X));
xlabel('time (centiseconds)');
ylabel('x(t)');
title(sprintf('x(t) vs. time, K = %d',K));

% Discussion
% You can see as K increases the corner of the boxcar waves become more and
% more well defined, and where K is close to 1, is behaves more
% sinusoidally. This is because as more and more sinusoids are combined to
% create the function the edges become more defined, but you can still see
% the large oscillations that happen near the edge with some of the plots
% (i.e. k = 100, k = 10)

type('pulseTrainDFS.m');
type('sumFourierSeries.m');
function plotComplexData(t,fterm)
% function plotComplexData(t,fterm)
% This function takes as an input a time vector and a calculated fourier
% vector, and plots them. It creates 4 plots:
%   - the real part of fterm vs. time
%   - the imaginary part of fterm vs. time
%   - absolute value vs. time
%   - unwrapped phase

% Plot the real vs time
subplot(2,2,1);
plot(t,real(fterm));
title('Real(fterm) vs. time');
ylabel('Real(fterm)');
xlabel('Time (sec)');
% Plot the imaginary vs time
subplot(2,2,2);
plot(t,imag(fterm));
title('Imaginary(fterm) vs. time');
ylabel('Imaginary(fterm)');
xlabel('Time (sec)');
% Plot the absolute value vs time
subplot(2,2,3);
plot(t,abs(fterm));
title('AbsoluteValue(fterm) vs. time');
ylabel('AbsoluteValue(fterm)');
xlabel('Time (sec)');
% Plot the unwrapped phase
subplot(2,2,4);
plot(unwrap(angle(fterm)));
title('Unwrapped Phase vs time');
ylabel('Unwrapped phase (Rad)');
xlabel('Time (centi-sec)');
suptitle('Complex Data Plots for fTerm');
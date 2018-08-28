function C = myCWT(data,fs,psi,scaleVec)
% function C = myCWT(data,fs,psi,scaleVec)
% Inputs:
%    'data' is vector of input data, length numSamps
%    fs is sampling rate of 'data' vector
%    psi is precomputed matrix of wavelet shapes, size
%        L x numScales, where L is the length of each wavelet
%    scaleVec is vector of scales used to compute psi, length numScales
% Output:
%    C:  numSamps x numScales CWT output

numScales = length(scaleVec);
% Ask Prof. Tracey about this in office hours
d = size(psi);
L = d(1);

% file loops over scales.  At each scale, convolves 'data' with a
% the apprpriate column of psi, and stores result in correspoding
% column of C

for ss = 1:numScales
    h = conj(psi(:,ss))/(fs*scaleVec(ss));
    C(:,ss) = conv(data,h,'same');
end

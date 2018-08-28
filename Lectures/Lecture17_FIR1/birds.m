
bLoon = false;
if bLoon
	[snd,fs]=audioread('loonwail.wav');
else
	[snd,fs]=audioread('chicken.wav');
	snd = snd(:,1);
end

nfft = 2048;
noverlap = round(nfft*.95);
%win = boxcar(nfft);
win = blackman(nfft);
%win = hanning(nfft);
figure
spectrogram(snd,win,noverlap,nfft,fs);

if bLoon,
	caxis([-80 -10])
	title('loon call')
	%xlim([0 10000])
	
else
	caxis([-90 -10])
	xlim([0 3])
	title('chicken clucking')	

end

%soundsc(snd,fs)

colorbar

disp('Loading data')
load('/Users/emilyschlafly/BU/DATA/BW09/BW09_Seizure1.mat')
mea = Neuroport;
ecog = ECoG;
clear Neuroport ECoG

% Preprocessing
tmax = Inf;

% Fs = ECoG.SamplingRate;
% interval = Neuroport.SamplingRate/ECoG.SamplingRate;
% time = ECoG.Time;
% tmax = min([tmax, length(ECoG.Time), floor(length(Neuroport.Time)/interval)]);

%% Filter
bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
	'CutoffFrequency1',2,'CutoffFrequency2',50, ...
	'SampleRate', mea.SamplingRate);
mea.lfp = single(filtfilt(bpFilt, double(mea.Data)));
mea.lfp(:, mea.BadChannels) = [];

bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
	'CutoffFrequency1',2,'CutoffFrequency2',50, ...
	'SampleRate', ecog.SamplingRate);
ecog.lfp = single(filtfilt(bpFilt, double(ecog.Data)));
ecog.lfp(:, ecog.BadChannels) = [];

bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
	'CutoffFrequency1',3e2,'CutoffFrequency2',3e3, ...
	'SampleRate',mea.SamplingRate);
mea.mua = single(filtfilt(bpFilt, double(mea.Data)));
mea.mua(:, mea.BadChannels) = [];

bpFilt = designfilt('bandpassfir', 'FilterOrder', 150, ...
	'CutoffFrequency1', 50, 'CutoffFrequency2', 300, ...
	'SampleRate', mea.SamplingRate);
mea.highg = single(filtfilt(bpFilt, double(mea.Data)));
mea.highg(:, mea.BadChannels) = [];

bpFilt = designfilt('bandpassfir', 'FilterOrder', 150, ...
	'CutoffFrequency1', 50, 'CutoffFrequency2', 124, ...  % 150 in the Schevon paper, but max is 124.93 given sampling freq
	'SampleRate', ecog.SamplingRate);
ecog.highg = single(filtfilt(bpFilt, double(mea.Data)));
ecog.highg(:, ecog.BadChannels) = [];

mea.X = mea.ElectrodeXY(:, 1);
mea.X(mea.BadChannels) = [];

mea.Y = mea.ElectrodeXY(:, 2);
mea.Y(mea.BadChannels) = [];

ecog.X = ecog.ElectrodeXY(:, 1);
ecog.X(ecog.BadChannels) = [];
ecog.Y = ecog.ElectrodeXY(:, 2);
ecog.Y(ecog.BadChannels) = [];

%% Get MUA events
% ?Negative peaks in this signal were detected and those peaks that
% exceeded four times the s.d. of the signal in the negative direction, and
% that occurred >1ms after immediately preceding peaks, were retained as
% multiunit timestamps and waveforms (Smith et al., 2016)

intervalM = mea.SamplingRate / 1e3;  % samples per ms
mea.artefacts = abs(zscore(mea.mua)) > 8;
mea.mua(mea.artefacts) = 0;
temp = mea.mua;
temp(temp > 0) = 0;
mea.events = false(size(temp));
for ch = 1:size(temp, 2)
	temp = mea.mua(:, ch);
	[~, inds] = findpeaks(-zscore(temp), ...
		'minpeakdistance', intervalM, 'minpeakheight', 2);  % Few peaks if sd=4
	mea.events(inds, ch) = true;
end

%% Compute firing rate
% Firing rate was calculated from these multiunit timestamps in 100-ms
% windows every 25 ms (Smith et al., 2016)

window = 100 * intervalM;  % 100 ms
[T, N] = size(mea.mua);
ts = 1:(25*intervalM):(T-window);
firingRate = zeros(numel(ts), N, 'single');

for i = 1:length(ts)
	t = ts(i);
	inds = t:(t + window - 1);
	firingRate(i, :) = sum(mea.events(inds, :)) * 10;  % events per second
end

%% Compute high-gamma instantaneous amplitude
% High-g instantaneous amplitude was defined as the absolute value of the
% Hilbert transform of this signal. (Smith et al., 2016)

mea.hga = abs(hilbert(mea.highg));
ecog.hga = abs(hilbert(ecog.highg));


%% Compute LFP spectral content
FsM = mea.SamplingRate;
intervalM = round(FsM);
overlapM = round(FsM * .95);
nfftM = intervalM;

FsE = ecog.SamplingRate;
intervalE = round(FsE);
overlapE = round(FsE * .95);
nfftE = intervalE;

Slfp = struct();
Smua = struct();
Secog = struct();

[Slfp.S, Slfp.F, Slfp.T, Slfp.P] = spectrogram(mea.lfp(:, 1) - mean(mea.lfp(:, 1)), ...
	intervalM, overlapM, nfftM, FsM); 
[Smua.S, Smua.F, Smua.T, Smua.P] = spectrogram(mea.mua(:, 1) - mean(mea.mua(:, 1)), ...
	intervalM, overlapM, nfftM, FsM); 
[Secog.S, Secog.F, Secog.T, Secog.P] = spectrogram(ecog.Data(:, 1) - mean(ecog.Data(:, 1)), ...
	intervalE, overlapE, nfftE, FsE); 

%% Morlet transformation

inds = logical((mea.Time > 0) .* (mea.Time < 60));
[cfsM, fM] = cwt(double(mean(mea.lfp(inds, :), 2)), FsM, 'amor');
figure(2); imagesc(mea.Time(inds), fM(fM <= 50), 10 * log10(abs(cfsM(f <= 50, :))))
axis xy
colorbar

inds = logical((ecog.Time > -3) .* (ecog.Time < 63));
[cfs, f] = cwt(double(mean(ecog.lfp(inds, :), 2)), FsE, 'amor');
figure(3); imagesc(ecog.Time(inds), f(f <= 50), (conv2(abs(cfs(f <= 50, :)), ones(10, 1000))))
colormap('jet')
axis xy
colorbar

%% Plot spectrograms
figure(1); clf
colormap('jet')
subplot(221)
imagesc(Slfp.T, Slfp.F(Slfp.F < 50), 10 * log10(Slfp.P(Slfp.F <= 50, :)))
colorbar
axis xy
ylim([2 50])
ylabel('Frequency [Hz]')
title('LFP')


subplot(222)
imagesc(Smua.T, Smua.F, 10 * log10(Smua.P), [-40 0])
colorbar
ylim([300 3e3])
axis xy

title('MUA')
% ylabel('Frequency [Hz]')

subplot(223)
imagesc(Secog.T, Secog.F, 10 * log10(Secog.P))
colorbar
ylim([2 50])
axis xy
title('ECoG')

print(6, '6', '-dpng')

%% Define epochs
% The ictal wavefront epoch was defined as the time of the first channel?s
% mean minus its s.d. until the last channel?s mean plus its s.d. (Smith et
% al., 2016).

frSm = smoothdata(firingRate, 'gaussian', 100);
peakRate = zeros(size(firingRate, 2), 1);
frT = mea.Time(ts);
inds = logical((frT >= 0) .* (frT <= 20));
frT = frT(inds);
mea.epochs = zeros(3, 1);

for i = 1:numel(peakRate)
	[~, ind] = max((frSm(inds, i)));
	peakRate(i) = frT(ind);
end

peakRate([68 72]) = [];  % These two were not recruited... one is definitely some kind of inhibitory cell
[firstT, firstCh] = min(peakRate);
[lastT, lastCh] = max(peakRate);

firstT = firstT - std(peakRate);
lastT = lastT + std(peakRate);

mea.epochs(1) = firstT;  % Recruitment start time
mea.epochs(2) = lastT;  % Post recruitment start time

[~, mea.dischargeTimes] = findpeaks(-mean(mea.highg, 2), ...
	'minpeakheight', 500, 'minpeakprominence', 100);
mea.dischargeTimes = mea.Time(mea.dischargeTimes);

mea.dischargeTimes(mea.dischargeTimes < mea.epochs(2)) = [];  % exclude recruitment discharges
mea.dischargeTimes(mea.dischargeTimes > mea.EndTime) = [];  % exclude post termination discharges

mea.IDIs = diff(mea.dischargeTimes);
mea.cov = zeros(floor((numel(mea.IDIs) - 10) / 3), 1);
for i = 1:numel(mea.cov)
	ii = 3*(i-1);
	interval = min(30, numel(mea.IDIs) - ii);
	temp = mea.IDIs(ii+1:ii+interval);
	mea.cov(i) = std(temp) / mean(temp);
end

[~, mea.epochs(3)] = max(diff(mea.cov));
mea.epochs(3) = mea.dischargeTimes(3 * (mea.epochs(3) - 1) + 1);
figure(10); plot(diff(mea.cov))


%%
dataNPds.lfp = interp1(Neuroport.Time, mea.lfp, time, 'nearest', 'extrap');
dataNPds.mua = interp1(Neuroport.Time, mea.mua, time, 'nearest', 'extrap');

dataE = ECoG.Data;
dataE(:, ECoG.BadChannels) = [];

xyE = ECoG.ElectrodeXY;
xyE(ECoG.BadChannels, :) = [];

xyNP = Neuroport.ElectrodeXY;
xyNP(Neuroport.BadChannels, :) = [];

xE = xyE(:, 1);
yE = xyE(:, 2);

xNP = xyNP(:, 1);
yNP = xyNP(:, 2);

chE = size(dataE, 2);
chNP = size(mea, 2);


%% Whiten
disp('Whitening')
[~, startInd] = min(abs(time));
[~, endInd] = min(abs(time - (time(end) - ECoG.Padding(2))));
dataNPW.mua = [zca_whitening(mea.mua(1:startInd-1, :)); ...
	zca_whitening(mea.mua(startInd:endInd, :)); ...
	zca_whitening(mea.mua(endInd+1:end, :))];
dataNPW.lfp = [zca_whitening(mea.lfp(1:startInd-1, :)); ...
	zca_whitening(mea.lfp(startInd:endInd, :)); ...
	zca_whitening(mea.lfp(endInd+1:end, :))];
dataEW = [zca_whitening(dataE(1:startInd-1, :)); ...
	zca_whitening(dataE(startInd:endInd, :)); ...
	zca_whitening(dataE(endInd+1:end, :))];

%% Remove 60 Hz signal
bpFilt = designfilt('bandpassfir','FilterOrder',20, ...
	'CutoffFrequency1',59.5,'CutoffFrequency2',60.5, ...
	'SampleRate',FsM);
dataEWf = dataEW - filtfilt(bpFilt, dataEW);
dataNPWf = dataNPW - filtfilt(bpFilt, double(dataNPW));

%% Reshape
disp('Reshaping')

dataEWr = nan(tmax, max(xE), max(yE));
dataNPWr = nan(tmax, max(xNP), max(yNP));

sE = size(dataEWr);
sNP = size(dataNPWr);


for t = 1:tmax
	dataEWr(t + (xE - 1) * sE(1) + (yE - 1) * prod(sE(1:2))) = ...
		dataEW(t, :);
	dataNPWr(t + (xNP - 1) * sNP(1) + (yNP - 1) * prod(sNP(1:2))) = ...
		dataNPW(t, :);

end


%%
dataNPWsm = reshape(smooth(dataNPW), [], size(dataNPW, 2)); %dataNP = dataNP - min(dataNP(:)) + 1;
dataEWsm = reshape(smooth(dataEW), [], size(dataEW, 2)); %dataE = dataE - min(dataE(:)) + 1;
%% Create a video of the whitened data

if 0
	disp('Making video of whitened data')
	CREATEVID = false;
	if CREATEVID
		v = VideoWriter('vis', 'MPEG-4');
		v.FrameRate = Fs;
		open(v);
	end

	mNP = quantile(dataNPWsm(:), .05);
	mE = quantile(dataEWsm(:), .05);

	MNP = quantile(dataNPWsm(:), .95);
	ME = quantile(dataEWsm(:), .95);

	map = make_diverging_colormap('cool', 1);

	try
		close 1;
	catch MException
	end
	for t = 1:1:tmax
		if t < startInd
			desc = 'preictal';
		elseif t > endInd
			desc = 'postictal';
		else
			desc = 'ictal';
		end
		figure(1);
		set(1, 'position', [358   443   778   254]);
		colormap(map)
		p1 = subplot(1,5,1:3);
		scatter(xE, yE, 150, dataEWsm(t, :), 's', 'filled')
		set(gca, 'clim', [mE ME], 'Color', .15*[1 1 1]);
		axis image
		xlim([0 max(xE)+1])
		ylim([0 max(yE)+1])
	% 	imagesc(squeeze(dataE(t, :, :)), [mE ME]);
		colorbar
		title(sprintf('T = %0.2f (%s)', ECoG.Time(t), desc))


		p2 = subplot(1, 5, 4:5);
		scatter(xNP, yNP, 200, dataNPWsm(t, :), 's', 'filled')
		set(gca, 'clim', [mNP MNP], 'Color', .15*[1 1 1]);
		axis image
		xlim([0 max(xNP)+1])
		ylim([0 max(yNP)+1])
	% 	imagesc(squeeze(dataNP(t, :, :)), [mNP MNP]);
		colorbar


		if CREATEVID
			frame = getframe(gcf);
			writeVideo(v,frame);
		else
			drawnow()
		end

		if ~mod(t, 10)
			disp(['t = ', num2str(t)])
		end
	end

	if CREATEVID
		close(v)
		disp('Video saved and closed')
	end
end
%% Create a video of the vector field of the whitened data

if 1
	disp('Making video of vector field')
	CREATEVID = true;
	if CREATEVID
		v = VideoWriter('vectors');
		v.FrameRate = FsM;
		open(v);
	end

	mNP = quantile(dataNPWsm(:), .05);
	mE = quantile(dataEWsm(:), .05);

	MNP = quantile(dataNPWsm(:), .95);
	ME = quantile(dataEWsm(:), .95);

	map = make_diverging_colormap('cool', 1);

	xxE = 1:max(xE); 
	yyE = 1:max(yE);

	xxNP = 1:max(xNP);
	yyNP = 1:max(yNP);

	try
		close 9;
	catch MException
	end

	convMat = -1/8 * ones(3); convMat(2, 2) = 1;
	for t = 1:1:tmax
		if t < startInd
			desc = 'preictal';
		elseif t > endInd
			desc = 'postictal';
		else
			desc = 'ictal';
		end
		figure(9);
		set(9, 'position', [358   443   778   254]);
		colormap(map)
		p1 = subplot(1,5,1:3);
		temp = squeeze(dataEWr(t, :, :));
		temp(isnan(temp)) = 0; temp = conv2(temp, convMat, 'same');
		contour(xxE, yyE, temp', linspace(mE, ME, 10)); hold on
		[fx, fy] = gradient(temp);
		quiver(xxE, yyE, fx', fy', 'LineWidth', 2, 'color', [1 1 1]); hold off;
		set(gca, 'clim', [mE ME], 'Color', .15*[1 1 1]);
		axis image
		xlim([0 max(xE)+1])
		ylim([0 max(yE)+1])
	% 	imagesc(squeeze(dataE(t, :, :)), [mE ME]);
		colorbar
		title(sprintf('T = %0.2f (%s)', ECoG.Time(t), desc))


		p2 = subplot(1, 5, 4:5);
		temp = squeeze(dataNPWr(t, :, :));
		temp(isnan(temp)) = 0; temp = conv2(temp, convMat, 'same');
		[fx, fy] = gradient(temp);
		contour(xxNP, yyNP, temp, linspace(mNP, MNP, 10)); hold on
		quiver(xxNP, yyNP, fx, fy, 'LineWidth', 2, 'color', [1 1 1]); hold off;	
		set(gca, 'clim', [mNP MNP], 'Color', .15*[1 1 1]);
		axis image
		xlim([0 max(xNP)+1])
		ylim([0 max(yNP)+1])
	% 	imagesc(squeeze(dataNP(t, :, :)), [mNP MNP]);
		colorbar


		if CREATEVID
			frame = getframe(gcf);
			writeVideo(v,frame);
		else
			drawnow()
		end

		if ~mod(t, 10)
			disp(['t = ', num2str(t)])
		end
	end

	if CREATEVID
		close(v)
		disp('Video saved and closed')
	end
end

%% Create a LIC video

if input('Make LIC video? (0/1) ')
	lic;
end

%% Choose which channels and times to plot
rows = 1:tmax;  % times
cols = [53 61];  % channels
%% Plot raw data
figure(2); fullwidth(true)
subplot(211)
plot(zscore((mea(rows, cols))) + (1:numel(cols)) * 5);
axis('tight')
title('Neuroport (zscored)')

subplot(212)
plot(zscore((dataE(rows, cols))) + (1:numel(cols)) * 5);
axis('tight')
title('ECoG (zscored)')

%% Plot whitened data
figure(3); fullwidth(true)
p1 = subplot(211);
plot(time(rows), medfilt1(dataNPW(rows, cols), 10) + (1:numel(cols)) * 2*max(var(dataNPW(rows, cols))));
axis('tight')
title('Neuroport (whitened and smoothed)')

p2 = subplot(212);
plot(time(rows), medfilt1(dataEW(rows, cols), 10) + (1:numel(cols)) * 2*max(var(dataEW(rows, cols))));
axis('tight')
title('ECoG (whitened and smoothed)')
linkaxes([p1 p2], 'x')
%% Plot difference between whitened and raw data
figure(4); fullwidth(true)
subplot(211)
plot(zscore(mea(rows, cols)) - zscore(dataNPW(rows, cols)) + ...
	(1:numel(cols)) * 5)

subplot(212)
plot(zscore(dataE(rows, cols)) - zscore(dataEW(rows, cols)) + ...
	(1:numel(cols)) * 5)

%% Correlation
NLAGS = 1000;
nChanE = size(dataE, 2);
nChanNP = size(mea, 2);
ccfE = reshape(xcorr(dataEWf, NLAGS, 'coeff'), [], nChanE, nChanE);
ccfNP = reshape(xcorr(dataNPWf, NLAGS, 'coeff'), [], nChanNP, nChanNP);

%% Spectrogram (compute)
% Se.Fs = ECoG.SamplingRate;
% Snp.Fs = ECoG.SamplingRate;  % remember that you downsampled above

intervalM = round(Se.Fs);
overlapM = round(Se.Fs * .95);
nfftM = intervalM;
[~, Se.F, Se.T, ~] = spectrogram(dataEWf(:, 1), ...
			intervalM, overlapM, nfftM, FsM); 
Se.P = zeros(numel(Se.F), numel(Se.T), chE);
Snp.P = zeros(numel(Se.F), numel(Se.T), chNP);

for ch=1:chNP
	if 1e6 <= chE
		[Se.S, Se.F, Se.T, Se.P(:, :, ch)] = spectrogram(dataEW(:, ch), ...
			intervalM, overlapM, nfftM, FsM); 
	end

	[Snp.S, Snp.F, Snp.T, Snp.P(:, :, ch)] = spectrogram(mea.mua(:, ch), ...
		intervalM, overlapM, nfftM, FsM); 
end


%% Spectrogram (plot)

% for ch=1:chNP
figure(6); fullwidth(); clf
colormap('parula')
for ch = [2 3 2 3]
	
	if 1e5 <= chE
		subplot(211)
		imagesc(Se.T - 60, Se.F, 10 * (log10(conv2(Se.P(:, :, ch), ones(3), 'same'))))
		set(gca, 'clim', [-60 0]);
		hold on; plot(time(ECoG.SyncInfo.EventIdx), Se.F(61) * ones(6, 1), 'r*'); hold off
		colorbar
		axis xy
		ylabel('Frequency [Hz]')
		title(['ECoG ch' num2str(ch)])
	end

	subplot(212)
	imagesc(Snp.T - 60, Snp.F, 10 * log10(conv2(Snp.P(:, :, ch), ones(3))))
% 	set(gca, 'clim', [-60 0]);
	colorbar, axis xy, ylabel('Frequency [Hz]');
	title(['Neuroport ch' num2str(ch)])

	drawnow()
	pause(1)
end

%% Spectrogram (summary)
figure(61)
subplot(211)
imagesc(Se.T, Se.F, 10 * log10(var(Se.P, [], 3)))
set(gca, 'clim', [-60 0]); colorbar, axis('xy'), ylabel('Frequency [Hz]');
% ylim([0 20])
title('ECoG')

subplot(212)
imagesc(Snp.T, Snp.F, 10 * log10(var(Snp.P, [], 3)))
set(gca, 'clim', [-60 0]); colorbar, axis('xy'), ylabel('Frequency [Hz]');
% ylim([0 20])
title('Neuroport')



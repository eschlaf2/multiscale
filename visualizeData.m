
disp('Loading data')
% load('/Users/emilyschlafly/BU/Work/DATA/BW09/BW09_Seizure1.mat')

tmax = Inf;

Fs = ECoG.SamplingRate;
interval = Neuroport.SamplingRate/ECoG.SamplingRate;
time = ECoG.Time;
tmax = min([tmax, length(ECoG.Time), floor(length(Neuroport.Time)/interval)]);

dataNP = interp1(Neuroport.Time, single(Neuroport.Data), time, 'nearest', 'extrap');
dataNP(:, Neuroport.BadChannels) = [];

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
chNP = size(dataNP, 2);


%% Whiten
disp('Whitening')
[~, startInd] = min(abs(time));
[~, endInd] = min(abs(time - (time(end) - ECoG.Padding(2))));
dataNPW = [zca_whitening(dataNP(1:startInd-1, :)); ...
	zca_whitening(dataNP(startInd:endInd, :)); ...
	zca_whitening(dataNP(endInd+1:end, :))];
dataEW = [zca_whitening(dataE(1:startInd-1, :)); ...
	zca_whitening(dataE(startInd:endInd, :)); ...
	zca_whitening(dataE(endInd+1:end, :))];

%% Remove 60 Hz signal
bpFilt = designfilt('bandpassfir','FilterOrder',20, ...
	'CutoffFrequency1',59.5,'CutoffFrequency2',60.5, ...
	'SampleRate',ECoG.SamplingRate);
dataEW(:, 1:4) = dataEW(:, 1:4) - filtfilt(bpFilt, dataEW(:, 1:4));

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
		v.FrameRate = Fs;
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
rows = 1:3.5e4;  % times
cols = [1 47:49];  % channels
%% Plot raw data
figure(2); fullwidth(true)
subplot(211)
plot(zscore((dataNP(rows, cols))) + (1:numel(cols)) * 5);
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
plot(zscore(dataNP(rows, cols)) - zscore(dataNPW(rows, cols)) + ...
	(1:numel(cols)) * 5)

subplot(212)
plot(zscore(dataE(rows, cols)) - zscore(dataEW(rows, cols)) + ...
	(1:numel(cols)) * 5)

%% Correlation
nChanE = size(dataE, 2);
nChanNP = size(dataNP, 2);
ccfE = reshape(xcorr(dataEW, 300, 'coeff'), [], nChanE, nChanE);
ccfNP = reshape(xcorr(dataNPW, 300, 'coeff'), [], nChanNP, nChanNP);

%% Spectrogram (compute)
% Se.Fs = ECoG.SamplingRate;
% Snp.Fs = ECoG.SamplingRate;  % remember that you downsampled above

interval = round(Se.Fs);
overlap = round(Se.Fs * .95);
nfft = interval;
[~, Se.F, Se.T, ~] = spectrogram(dataEW(:, 1), ...
			interval, overlap, nfft, Fs); 
Se.P = zeros(numel(Se.F), numel(Se.T), chE);
Snp.P = zeros(numel(Se.F), numel(Se.T), chNP);

for ch=1:chNP
	if ch <= chE
		[Se.S, Se.F, Se.T, Se.P(:, :, ch)] = spectrogram(dataEW(:, ch), ...
			interval, overlap, nfft, Fs); 
	end

	[Snp.S, Snp.F, Snp.T, Snp.P(:, :, ch)] = spectrogram(dataNPW(:, ch), ...
		interval, overlap, nfft, Fs); 
end


%% Spectrogram (plot)

for ch=1:chNP
	if ch <= chE
		figure(6); fullwidth(); clf
		colormap('parula')
		subplot(211)
		imagesc(Se.T, Se.F, 10 * log10(conv2(Se.P(:, :, ch), ones(3))))
		set(gca, 'clim', [-60 0]);
		colorbar
		axis xy
		ylabel('Frequency [Hz]')
		title('ECoG')
	end

	subplot(212)
	imagesc(Snp.T, Snp.F, 10 * log10(conv2(Snp.P(:, :, ch), ones(3))))
	set(gca, 'clim', [-60 0]), colorbar, axis xy, ylabel('Frequency [Hz]'), title('Neuroport')

	drawnow()
	pause(1e-2)
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



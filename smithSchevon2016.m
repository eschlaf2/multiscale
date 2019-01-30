disp('Loading data')
load('/Users/emilyschlafly/BU/DATA/BW09/BW09_Seizure1.mat')
if ~exist('mea', 'var')
	mea = Neuroport;
	ecog = ECoG;
	clear Neuroport ECoG
end

set(groot, 'defaultLineLineWidth', 1);
set(groot, 'defaultAxesLineWidth', 1);
set(groot, 'defaultFigureColor', [1 1 1], 'defaultAxesFontSize', 16, 'defaultHistogramLineWidth',1);
% set(groot, 'defaultAxesFontWeight', 'bold');
set(groot, 'defaultAxesFontSize', 16);
c = lines(7);
aic = @(predictors, dev) dev + 2 * size(predictors, 2);

tmax = Inf;

%% Filter
disp('Filtering');
if ~isfield(mea, 'lfp')

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
	ecog.highg = single(filtfilt(bpFilt, double(ecog.Data)));
	ecog.highg(:, ecog.BadChannels) = [];

	mea.X = mea.ElectrodeXY(:, 1);
	mea.X(mea.BadChannels) = [];

	mea.Y = mea.ElectrodeXY(:, 2);
	mea.Y(mea.BadChannels) = [];

	ecog.X = ecog.ElectrodeXY(:, 1);
	ecog.X(ecog.BadChannels) = [];
	ecog.Y = ecog.ElectrodeXY(:, 2);
	ecog.Y(ecog.BadChannels) = [];
end


%% Whiten data
% mea.muaW = zca_whitening(mea.mua);
% mea.lfpW = zca_whitening(mea.lfp);
% mea.highgW = zca_whitening(mea.highg);
% ecog.highgW = zca_whitening(ecog.highg);
% ecog.lfpW = zca_whitening(ecog.lfp);

%% Get MUA events
% Negative peaks in this signal were detected and those peaks that
% exceeded four times the s.d. of the signal in the negative direction, and
% that occurred >1ms after immediately preceding peaks, were retained as
% multiunit timestamps and waveforms (Smith et al., 2016)

disp('Getting MUA events');

if ~isfield(mea, 'events')
	mea = mua_events(mea);
end

%% Compute firing rate
% Firing rate was calculated from these multiunit timestamps in 100-ms
% windows every 25 ms (Smith et al., 2016)

disp('Computing firing rate');
if ~isfield(mea, 'firingRate')
	mea = mua_firing_rate(mea);
end

%% Order traces by firing rate

order = zeros(size(firingRate, 2), 1);
tt = 2400:3000;
for i = 1:numel(order)
	[~, order(i)] = max((firingRate(tt, i)));
end

[~, ord] = sort(order);

%% Compute high-gamma instantaneous amplitude
% High-g instantaneous amplitude was defined as the absolute value of the
% Hilbert transform of this signal. (Smith et al., 2016)

disp('Computing high-g');
if ~isfield(mea, 'hga')
	mea.hga = abs(hilbert(mea.highgW));
	ecog.hga = abs(hilbert(ecog.highgW));
end

%% Compute LFP spectral content

disp('Computing spectral content')
if ~exist('Slfp', 'var')
	FsM = mea.SamplingRate;
	samplingRateMs = round(FsM);
	overlapM = round(FsM * .95);
	nfftM = samplingRateMs;

	FsE = ecog.SamplingRate;
	intervalE = round(FsE);
	overlapE = round(FsE * .95);
	nfftE = intervalE;

	Slfp = struct();
	Smua = struct();
	Secog = struct();

	[Slfp.S, Slfp.F, Slfp.T, Slfp.P] = spectrogram(mea.lfpW(:, 1) - mean(mea.lfpW(:, 1)), ...
		samplingRateMs, overlapM, nfftM, FsM); 
	[Smua.S, Smua.F, Smua.T, Smua.P] = spectrogram(mea.muaW(:, 1) - mean(mea.muaW(:, 1)), ...
		samplingRateMs, overlapM, nfftM, FsM); 
	[Secog.S, Secog.F, Secog.T, Secog.P] = spectrogram(ecog.lfpW(:, 1) - mean(ecog.lfpW(:, 1)), ...
		intervalE, overlapE, nfftE, FsE); 
end

%% Morlet transformation

disp('Computing Morlet transformation')
if ~exist('cfsM', 'var')
	step = 10;
% 	indsM = logical((mea.Time > 0) .* (mea.Time < 60));
	[cfsM, fM] = cwt(double(mean(mea.lfpW(1:step:end, :), 2)), FsM, 'amor');
	figure(2); 
	subplot(1, 2, 1); imagesc(mea.Time(1:step:end), fM(fM <= 50), 10 * log10(abs(cfsM(fM <= 50, :))))
	title('MEA spectrum (dB)'); colormap('parula')
	axis xy
	xlabel('Time (s)')
	ylabel('Frequency (Hz)')
	colorbar

	inds = logical((ecog.Time > 0) .* (ecog.Time < 60));
% 	tt = ecog.Time(inds);
% 	inds2 = logical((tt > 0) .* (tt < 60));
	[cfs, f] = cwt(double(mean(ecog.lfpW(:, :), 2)), FsE, 'amor');
% 	figure(3); imagesc(ecog.Time(inds), f(f <= 50), (conv2(abs(cfs(f <= 50, :)), ones(10, 1000))))
	subplot(1, 2, 2);
	imagesc(ecog.Time(:), f(f <= 50), 10 * log10(abs(cfs(f <= 50, :))))
% 	colormap('jet')
	title('ECoG spectrum (dB)')
	xlabel('Time (s)')
	xticks(10:10:50)
	axis xy
	colorbar
end

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

%% Define epochs and IDIs
% The ictal wavefront epoch was defined as the time of the first channel?s
% mean minus its s.d. until the last channel?s mean plus its s.d. (Smith et
% al., 2016).

disp('Defining epochs')

firingRateSm = smoothdata(firingRate, 'gaussian', 100);
peakRate = zeros(size(firingRate, 2), 1);
frT = mea.Time(ts);
inds = logical((frT >= 0) .* (frT <= 20));  % Look in the first 20 seconds of the seizure
frT = frT(inds);
mea.epochs = zeros(3, 1);

for i = 1:numel(peakRate)
	[~, ind] = max((firingRateSm(inds, i)));
	peakRate(i) = frT(ind);
end

peakRate([68 72]) = [];  % These two electrodes look very different... one is definitely some kind of inhibitory cell
[firstT, firstCh] = min(peakRate);
[lastT, lastCh] = max(peakRate);

firstT = firstT - std(peakRate);
lastT = lastT + std(peakRate);

mea.epochs(1) = firstT;  % Recruitment start time
mea.epochs(2) = lastT;  % Post recruitment start time

% Discharge times: find peaks in firing rate
sd = std(mean(firingRate(mea.Time(ts) < 0, :), 2));
mn = mean(mean(firingRate(mea.Time(ts) < 0, :), 2));
[~, temp] = findpeaks((mean(firingRate, 2)) - mn, ...
	'minpeakheight', 2*sd, 'minpeakprominence', .25 * sd);
hga = (mean(mea.hga, 2) - mean(mean(mea.hga(mea.Time < 0, :), 2))) / ...
	std(mean(mea.hga(mea.Time < 0, :), 2));
hga = hga(ts(temp));  % high-gamma at discharge times
mea.dischargeInds = ts(temp(hga > 2));
mea.dischargeTimes = mea.Time(mea.dischargeInds);  % peaks with high-gamma at least one sd above the mean

mea.dischargeTimes(mea.dischargeTimes < mea.epochs(2)) = [];  % exclude recruitment discharges
mea.dischargeTimes(mea.dischargeTimes > (mea.Time(end) - mea.Padding(2))) = [];  % exclude post termination discharges

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
mea.dischargeTimes = mea.Time(mea.dischargeInds);
% figure(10); plot(diff(mea.cov))

%% Traveling wave direction

disp('Computing travelling wave direction')

ecog.dischargeInds = interp1(ecog.Time, 1:length(ecog.Time), mea.dischargeTimes, 'nearest');

Zm = wave_dir(mea);
Ze = wave_dir(ecog);

%%
figure(20); 
inds = logical((ecog.Time > -10) .* (ecog.Time < te));

Zma = interp1(mea.dischargeTimes, abs(Zm), ecog.Time(inds), 'linear'); 
Zmphi = interp1(mea.dischargeTimes, angle(Zm), ecog.Time(inds), 'linear');
subplot(2, 1, 1); [hAx,hLine1,hLine2] = plotyy(ecog.Time(inds), smooth(Zma, 100), ecog.Time(inds), smooth(Zmphi, 100)); axis(hAx, 'tight');
% hLine1.LineStyle = 'none'; hLine2.LineStyle = 'none';
% hLine1.Marker = 'o'; hLine2.Marker = 'o';
title('MEA waves')
ylabel(hAx(1), 'Amplitude'); ylabel(hAx(2), 'Angle')

Zea = interp1(mea.dischargeTimes, abs(Ze), ecog.Time(inds), 'linear');
Zephi = interp1(mea.dischargeTimes, angle(Ze), ecog.Time(inds), 'linear');
subplot(2, 1, 2); [hAx,hLine1,hLine2] = plotyy(ecog.Time(inds), smooth(Zea, 100), ecog.Time(inds), smooth(Zephi, 100)); axis(hAx, 'tight');
% hLine1.LineStyle = 'none'; hLine2.LineStyle = 'none';
% hLine1.Marker = 'o'; hLine2.Marker = 'o';
title('ECoG waves')
ylabel(hAx(1), 'Amplitude'); ylabel(hAx(2), 'Angle')
xlabel('Time (s)')

% subplot(3, 1, 3); [hAx,hLine1,hLine2] = plotyy(ecog.Time(ecog.dischargeInds), abs(Ze) - abs(Zm), ecog.Time(ecog.dischargeInds), angle(Ze) - angle(Zm));
% hLine1.LineStyle = 'none'; hLine2.LineStyle = 'none';
% hLine1.Marker = 'o'; hLine2.Marker = 'o';
% title('Difference')
% ylabel(hAx(1), 'Amplitude'); ylabel(hAx(2), 'Angle')
% xlabel('Time (s)')

figure(21); 
subplot(121); rose(angle(Zm)); title('MEA')
subplot(122); rose(angle(Ze)); title('ECoG')
%% Coherence

[cohg, conf, fm, tm] = cohgramc_rand_pairs(mea);
[cohgE, confE, fE, tE] = cohgramc_rand_pairs(ecog);
% [~, ~, fm, tm] = cohgramc_rand_pairs(mea);

mnCohg = cohgE;
mnCohg(isnan(mnCohg)) = 0;
mnCohg = mean(mnCohg, 3);
cmap = parula(64);
% cmap = [1 1 1; cmap];
figure(10); 
subplot(3, 1, 1:2); 
ax1 = imagesc(t - 60, fm(fm > 2), 10 * log10(mnCohg(:, fm > 2))'); axis xy;
colormap(cmap)
ylabel('Freq (Hz)');
title('MEA Coherence (dB)')
xticklabels([])
colorbar('south');

subplot(3, 1, 3); plot(t - 60, 10 * log10(mean(mnCohg(:, fm > 2), 2)));
axis tight;
xlabel('Time (s)');

%% spike-field coherence
clear params

dt = 1/mea.SamplingRate;
TW = 20;
ntapers = 2*TW-1;
params.Fs = 1e3;
step = ceil(mea.SamplingRate / params.Fs);
params.Fs = mea.SamplingRate / step;
params.tapers = [TW, ntapers];
params.pad = -1;
params.trialave = 1;
params.fpass = [2 50];
params.err = [1 .005];

window = 10;  % seconds
winstep = 1;  % seconds


inds = mea.Time < (te + 10);
T = sum(inds);
inds = inds(1:step:floor(sum(inds)/step)*step);

temp = mea.events(1:floor(T/step)*step, :);
temp = reshape(temp, step, [], size(temp, 2));
temp = squeeze(sum(temp));

% [C2, ~, ~, Syy, Snn, f2] = coherencycpb(mea.hga(inds, 1), mea.events(inds, 1), params);
ch = size(mea.events, 2);
[C,~,~,~,~,t,f,~,confC,~]=cohgramcpb((mea.lfpW(inds, 1)),(temp(:, 1)),[window winstep],params);
cohgSF = zeros([size(C) ch]);
cohgSF(:, :, 1) = C; 
for i = 2:ch
	[C,~,~,~,~,~,~,~,confC,~]=cohgramcpb((mea.lfpW(inds, i)),(temp(:, i)),[window winstep],params);
	cohgSF(:, :, i) = C;
end

mnCohg = cohgSF;
mnCohg(isnan(mnCohg)) = 0;
mnCohg = mean(mnCohg, 3);

% mnCohg = C;
% cmap = parula(64);
% cmap = [1 1 1; cmap];
figure(10); clf
subplot(3, 1, 1:2); 
ax1 = imagesc(t - 60, f, 10 * log10(mnCohg)'); axis xy;
% colormap('parula')
ylabel('Freq (Hz)');
title('Spike-field Coherence (dB)')
xticklabels([])
colorbar('south');

subplot(3, 1, 3); plot(t - 60, 10 * log10(mean(mnCohg, 2)));
axis tight;
xlabel('Time (s)');

%% Correlations
ch = size(firingRate, 2);
ccmat = zeros(200, ch);
for i = 1:ch
	[cc, lags] = xcorr(firingRate(:, i), 200, 'coeff');
	ccmat(:, i) = cc(lags > 0);
end
figure(23); clf
subplot(7,2,[1 3 5 7]); imagesc(lags(lags > 0) * 25e-3, 1:ch, ccmat')
title('Autocorrelation coefficients')
ylabel('Channel')
colorbar('southoutside')
xticks([])
subplot(7, 2, [9 11]); plot(lags(lags > 0), mean(ccmat, 2))
xlabel('Time (s)')
axis tight

% figure(24); 
ccov = cov(firingRate);
ccov = triu(ccov); ccov = ccov - diag(diag(ccov));
subplot(7, 2, [4 6 8 10]); imagesc(ccov)
title('Covariance of firing rates')
colorbar


%% Phase-amplitude coupling

ch = 10;
[phi, phiInd] = sort(angle(hilbert(mea.lfp(:, ch))));
figure(1); subplot(211); plot(phi, mea.mua(phiInd, ch)); axis tight; title('MUA vs Ictal phase');
figure(1); subplot(212); plot(phi, mea.hga(phiInd, ch)); axis tight; title('High-\gamma vs Ictal phase');


%% Plot raster

temp = find(mea.Time >=0, 1):find(mea.Time >= 5);
tt = repmat(mea.Time(temp(:)), 1, N);
yy = repmat(1:N, numel(temp), 1);
figure(2); 
subplot(511); plot(mea.Time(temp), mean(mea.lfp(temp, :), 2)); axis tight;
xticks([]); yticks([]); ylabel('LFP');
subplot(512); plot(mea.Time(temp), mean(mea.highg(temp, :), 2)); axis tight;
xticks([]); yticks([]); ylabel('High-\gamma')
subplot(5, 1, 3:5); scatter(tt(mea.events(temp, :)), yy(mea.events(temp, :)), 6, 'filled');
yticks([]); ylabel('MUA'); xlabel('Time (s)')
axis tight;

%% Plot traces
te = mea.Time(end) - mea.Padding(2);
inds = logical((mea.Time < 70) .* (mea.Time > -10));
figure(99); clf; % fullwidth(true); 
subplot(311); plot(mea.Time(inds), zscore(mean(mea.lfp(inds, :), 2)));
axis tight
hold on;  % show epochs
lims = get(gca, 'ylim');
limx = get(gca, 'xlim');
patch([limx(1), 0, 0, limx(1)], [lims(1) lims(1) lims(2) lims(2)], [.5 .5 .5], 'facealpha', .3, 'edgecolor', 'none')
patch([te, limx(2), limx(2), te], [lims(1) lims(1) lims(2) lims(2)], [.5 .5 .5], 'facealpha', .3, 'edgecolor', 'none')
plot([mea.epochs mea.epochs]', repmat(lims, 3, 1)', ':', 'linewidth', 4); hold off;
title('Mean LFP'); 

subplot(312); 
plot(mea.Time(inds), zscore(mean(mea.hga(inds, :), 2)))
axis tight
hold on;  % show epochs
lims = get(gca, 'ylim');
limx = get(gca, 'xlim');
patch([limx(1), 0, 0, limx(1)], [lims(1) lims(1) lims(2) lims(2)], [.5 .5 .5], 'facealpha', .3, 'edgecolor', 'none')
patch([te, limx(2), limx(2), te], [lims(1) lims(1) lims(2) lims(2)], [.5 .5 .5], 'facealpha', .3, 'edgecolor', 'none')
plot([mea.epochs mea.epochs]', repmat(lims, 3, 1)', ':', 'linewidth', 4); hold off;
% , ...
% 	mea.dischargeTimes(mea.dischargeTimes < 70), ...
% 	zeros(1, sum(mea.dischargeTimes < 70)), 'r*'); 
title('Mean High-\gamma amplitude')
ylabel({'z-scored voltage (\mu V),'; 'amplitude, rate (events/sec)'; ''})

subplot(313); 
inds = logical((mea.Time(ts) < 70) .* (mea.Time(ts) > -10));
plot(mea.Time(ts(inds)), ...
	zscore(mean(firingRate(inds, :), 2))); axis tight
hold on;  % show epochs
lims = get(gca, 'ylim');
limx = get(gca, 'xlim');
patch([limx(1), 0, 0, limx(1)], [lims(1) lims(1) lims(2) lims(2)], [.5 .5 .5], 'facealpha', .3, 'edgecolor', 'none')
patch([te, limx(2), limx(2), te], [lims(1) lims(1) lims(2) lims(2)], [.5 .5 .5], 'facealpha', .3, 'edgecolor', 'none')
plot([mea.epochs mea.epochs]', repmat(lims, 3, 1)', ':', 'linewidth', 4); hold off;
title('Mean firing rate'); 
xlabel('Time (s)')

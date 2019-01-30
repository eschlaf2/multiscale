function [cohg, conf, f, t] = cohgramc_rand_pairs(mea)

TW = 20;
ntapers = 2*TW-1;
params.Fs = 1e3;
step = ceil(mea.SamplingRate / params.Fs);
params.Fs = mea.SamplingRate / step;
params.tapers = [TW, ntapers];
params.pad = -1;
params.trialave = 1;
params.fpass = [1 13];
params.err = [1 .005];

window = 10;  % seconds
winstep = 1;  % seconds

te = mea.Time(end) - mea.Padding(2);

inds = mea.Time < (te + 10);
inds = inds(1:step:end);

% [C2, ~, ~, S1, S2, f2] = coherencyc(mea.lfpW(inds, 1), mea.lfpW(inds, 2), params);
electrodes = randsample(length(mea.X), 20);
[C,phi,~,~,~,t,f,confC,phistd]=cohgramc((mea.lfpW(inds, 1)),(mea.lfpW(inds, 2)),[window winstep],params);
cohg = zeros([size(C), 100], 'single');
conf = zeros(100, 1);
end

function [] = internal()


for i = 1:10
	for j = 1:10
		e1 = electrodes(i); e2 = electrodes(j + 10);
		[C,phi,~,~,~,~,~,confC,~] = ...
			cohgramc( (mea.lfpW(inds, e1)), (mea.lfpW(inds, e2)),...
			[window winstep], params);
% 		temp = C;
		conf((i - 1) * 10 + j) = confC;
		cohg(:, :, (i - 1) * 10 + j) = C;
	end
end
end
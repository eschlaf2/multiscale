
p.mE = interp1(predictor.morletEt, predictor.morletE, mea.Time, 'linear', 0)';
p.mM = interp1(predictor.morletMt, predictor.morletM, mea.Time, 'linear', 0)';
p.cohM = interp1(predictor.cohMt, predictor.cohM, mea.Time, 'linear', 0)';
p.cohE = interp1(predictor.cohEt, predictor.cohE, mea.Time, 'linear', 0)';
p.cohSF = interp1(predictor.cohSFt, predictor.cohSF, mea.Time, 'linear', 0)';
p.history = interp1(predictor.historyT, predictor.history(:, 1), mea.Time, 'linear', 'extrap')';
p.history = repmat(p.history(:), 1, 10);
for i = 1:9
	p.history(:, i+1) = circshift(p.history(:, 1), round(i * 25e-3 * mea.SamplingRate));
end
p.ZrhoM = interp1(predictor.Zt, predictor.ZrhoM, mea.Time, 'linear', 0)';
p.ZrhoE = interp1(predictor.Zt, predictor.ZrhoE, mea.Time, 'linear', 0)';
p.ZphiM = interp1(predictor.Zt, predictor.ZphiM, mea.Time, 'linear', 0)';
p.ZphiE = interp1(predictor.Zt, predictor.ZphiE, mea.Time, 'linear', 0)';
p.epoch = zeros(length(mea.Time), 1);
p.epoch(mea.Time > 0) = 1;
p.epoch(mea.Time > mea.epochs(1)) = 2;
p.epoch(mea.Time > mea.epochs(2)) = 3;
p.epoch(mea.Time > mea.epochs(3)) = 4;
p.epoch(mea.Time > te) = 5;
p.epoch = p.epoch + 1;

predictors = struct2array(p);
%%
numCh = size(meaEvents, 2);
time = (1:length(predictors))' * 1/3e4 - 60;
predictors = [time p.mM p.mE p.cohSF p.history(:, [1 7 9 10])];


numP = size(predictors, 2);
betas = zeros(numP + 1, numCh, 5);
% b = zeros(numCh, 1);
dev = zeros(numCh, 5);
statp = Inf(size(betas));

%%
for i = 1:numCh  %  channels
	for j = 1:6  %  epochs
		[betas(:, i, j), dev(i, j), stats] = glmfit(predictors(p.epoch == j, :), meaEvents(p.epoch == j, i), 'poisson');
		statp(:, i, j) = stats.p;
	end
	disp(i)
end

%%

temp = .05 - statp; 
% temp(:, [69 78], :) = [];
mask = statp == 0;
temp(mask) = nan;


%%
numCh = size(temp, 2);
figure(102);
for i = 1:numP+1
	datatemp = squeeze(temp(i, :, :));
	datatemp([69 72], :) = nan;
	maskp = squeeze(statp(i, :, :) < .05);
	datatemp = datatemp.*maskp;
	subplot(2, 5, i); plot(datatemp + .07*(1:5), '.', 'MarkerSize', 10); axis tight
	yticks([])
end
pos = [0.7957    0.1108    0.1153    0.3203];
lgd = legend('pre-seizure', 'pre-recruitment', 'ictal wavefront', 'post-recruitment', 'pre-termination'); lgd.Position = pos;

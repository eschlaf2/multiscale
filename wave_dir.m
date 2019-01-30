function Z = wave_dir(mea)

%%

halfWin = 25;  % Time (ms) on either side of a discharge
halfIntM = round(mea.SamplingRate * 1e-3 * halfWin);  % Hz * 1s/1000ms * halfWin ms (Time converted to samples)
nx = min(mea.X):max(mea.X);
ny = min(mea.Y):max(mea.Y);
[XX, YY] = meshgrid(nx, ny);
cmap = jet(length(mea.dischargeInds));
T = length(mea.dischargeInds);

figure(11); fullwidth; clf
% compass(Z); 
Z = zeros(T, 1, 'single');

DI = mea.dischargeInds(:);  % move x ms earlier
DI = [DI(1) - halfIntM; DI; DI(end) + halfIntM]; 
for i = 2:T+1
	window = (DI(i) - min(round((DI(i) - DI(i-1))), halfIntM): ...
		DI(i) + round(min((DI(i+1) - DI(i)) / 4, halfIntM/4)));
% 	temp = mea.lfp(window, :);
	[~, minSlope] = min(diff(mea.lfp(window, :)));
	mn = mean(minSlope(:));
	
	% convert to matrix
	weights = full(sparse(mea.X, mea.Y, minSlope));
	nullWeights = (weights == 0);
	weights(~nullWeights) = weights(~nullWeights) - mn;
	
	% interpolate missing electrodes if there are not only a few
	if sum(nullWeights(:)) < .1 * numel(weights)
		weights(nullWeights) = nan;
		tempH = zeros(size(weights));
		tempV = tempH;
		for j = nx
			mask = ~isnan(weights(j, :));
			tempH(j, :) = interp1(ny(mask), weights(j, mask), ny, 'linear', 'extrap');
		end
		for j = ny
			mask = ~isnan(weights(:, j));
			tempV(:, j) = interp1(nx(mask), weights(mask, j), nx, 'linear', 'extrap');
		end
		weights = (tempH + tempV) / 2;
		nullWeights = false(size(nullWeights));
	end
	
	wAll = weights;
	weights = wAll.*(wAll > 0);
% 	weights = wAll > 0;
	weights = (weights) / sum((weights(:)));  % normalize

	% late
	X = sum(weights(:) .* XX(:));
	Y = sum(weights(:) .* YY(:));
	
	% early
% 	weights = (weights - max(weights(:)));
% 	weights(nullWeights) = 0;
	

	weights = -wAll.*(wAll < 0);
% 	weights = wAll < 0;
	weights = weights / sum(weights(:));
	
	x = sum(weights(:) .* XX(:));
	y = sum(weights(:) .* YY(:));
	
	Z(i-1) = complex(X - x, Y - y);
	figure(11); subplot(221); p = compass(Z(i-1)); hold on;
	p.Color = cmap(i-1, :);
	p.LineWidth = 2;
	
% 	subplot(122); p = feather(Z);
% 	p.LineWidth = 2;
% 	drawnow()
% 	angle(Z)
	
	figure(11); set(11,'DefaultAxesColorOrder',cool(90))
	subplot(222); scatter(mea.X, mea.Y, 200, minSlope > mean(minSlope), 's', 'filled'); 
	hold on; scatter([x X]', [y Y]', 250, [0 1]', 'filled'); hold off; colorbar
	subplot(224); plot(1:numel(window), mea.lfp(window, :), '-', minSlope(:), mea.lfp(window(minSlope(:)), end), 'r*');
% 	pause();
	drawnow();
end

figure(11);
subplot(221); hold off;
subplot(223); feather(Z); axis tight
subplot(224); rose(angle(Z));

% function covariance_over_time(data, window)

data = dataNPW;
window = round(Fs);
[T, CH] = size(data);

c = zeros(CH, CH, T-window);
for i = 1:(T-window)
	temp = data(i:i+window-1, :);
	c(:, :, i) = cov(temp - mean(temp));
% 	figure(99);
% 	imagesc(c(:, :, i), [1e5 3.2e5])
% 	colorbar
% 	drawnow()
end

%%
mask = tril(ones(size(c(:, :, 1)))) - diag(ones(CH, 1));
mask = mask > 0;
cT = zeros(T-window, sum(mask(:)));
for i = 1:T-window
	temp = reshape(c(:, :, i), [], 1);
	cT(i, :) = temp(mask(:));
end
% cT = cell2mat(arrayfun(@(i) c(mask, i), 1:2, 'uni', 0))';
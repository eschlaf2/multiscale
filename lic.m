
CREATEVID = true;
if CREATEVID
	v = VideoWriter('LIC', 'MPEG-4');
	v.FrameRate = Fs;
	open(v);
end

mNP = quantile(dataNPWsm(:), .05);
mE = quantile(dataEWsm(:), .05);

MNP = quantile(dataNPWsm(:), .95);
ME = quantile(dataEWsm(:), .95);

imDim = 500;
[Xq,Yq] = meshgrid(1:imDim, 1:imDim);
% Xq = Xq(:); Yq = Yq(:);
noise = @() double(rand(imDim) >= .5);
dt = .25;
% 	msize = 3;
steps = 20;
alpha = .2;
map = (gray(20));

try
	close 10;
catch MException
end

convMat = -1/8 * ones(3); convMat(2, 2) = 1;
% 	Xtemp = repmat(Xq(:), steps+1, 1); Ytemp = repmat(Yq(:), steps+1, 1);
img = noise();
% 	warpCoords = sub2ind(size(Xq), Xq(:) + Xq(:).*fx*dt, Yq(:) + Yq(:).*fy*dt);
for t = 1:1:tmax
	if t < startInd
		desc = 'preictal';
	elseif t > endInd
		desc = 'postictal';
	else
		desc = 'ictal';
	end
	figure(10);
% 		set(10, 'position', [358   443   778   254]);
	colormap(map)

% 		p1 = subplot(1,5,1:3);
	temp = squeeze(dataEWr(t, :, :));
	temp(isnan(temp)) = 0; temp = conv2(temp, convMat, 'same');
	mask = (temp == 0);
% 	contour(xxE, yyE, temp', linspace(mE, ME, 10)); hold on
	[fx, fy] = gradient(temp);
	fxi = interp2(linspace(1, imDim, size(fx, 2)), ...
		linspace(1, imDim, size(fx, 1)), ...
		fx, Xq, Yq); 
	fyi = interp2(linspace(1, imDim, size(fy, 2)), ...
		linspace(1, imDim, size(fy, 1)), ...
		fy, Xq, Yq);
	maski = logical(interp2(linspace(1, imDim, size(mask, 2)), ...
		linspace(1, imDim, size(mask, 1)), ...
		mask, Xq, Yq, 'nearest'));
% 	Xtemp = [Xtemp(1+imDim^2:end); Xq(:)]; 
% 	Ytemp = [Ytemp(1+imDim^2:end); Yq(:)];
% 		imgtemp = noise();
	warpX = fxi(:)*dt; warpY = fyi(:)*dt;

	img = (1-alpha)*img + alpha*noise();
	for i = 1:steps
		warpCoords = sub2ind(size(img), mod(round(Xq(:)-1 + i*warpX(:)), imDim)+1, mod(round(Yq(:)-1 + i*warpY(:)), imDim) + 1);
		img(warpCoords) = .5*(img(:) + img(warpCoords));
	end
	img = (img - min(img(:)));
	img = img / max(img(:));
	img(img > .95) = 1;
	img(img < .05) = 0;
	imgtemp = img;
	imgtemp(maski) = nan;
	imgtemp = imgtemp - mean(imgtemp(:), 'omitnan') + .5;
% 		imgtemp(warpX + warpY == 0) = nan;
% 		img(warpX + warpY == 0) = nan;
	imagesc(imgtemp', [.25 .75]); axis xy; axis image
% 		hold on;
% 		xlim([-1 imDim]);
% 		ylim([-1 imDim]);
% 		drawnow()
% 		hold off
% 		end
% 	quiver(xxE, yyE, fx', fy', 'LineWidth', 2, 'color', [1 1 1]); hold off;
% 	set(gca, 'clim', [mE ME], 'Color', .15*[1 1 1]);
% 	axis image
% 	xlim([0 max(xE)+1])
% 	ylim([0 max(yE)+1])
% 	imagesc(squeeze(dataE(t, :, :)), [mE ME]);
% 		colorbar
	title(sprintf('T = %0.2f (%s)', ECoG.Time(t), desc))


% 		p2 = subplot(1, 5, 4:5);
% 		temp = squeeze(dataNPWr(t, :, :));
% 		temp(isnan(temp)) = 0; temp = conv2(temp, convMat, 'same');
% 		[fx, fy] = gradient(temp);
% 		contour(xxNP, yyNP, temp, linspace(mNP, MNP, 10)); hold on
% 		quiver(xxNP, yyNP, fx, fy, 'LineWidth', 2, 'color', [1 1 1]); hold off;	
% 		set(gca, 'clim', [mNP MNP], 'Color', .15*[1 1 1]);
% 		axis image
% 		xlim([0 max(xNP)+1])
% 		ylim([0 max(yNP)+1])
% 	% 	imagesc(squeeze(dataNP(t, :, :)), [mNP MNP]);
% 		colorbar


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

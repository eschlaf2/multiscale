function map = make_diverging_colormap(cmap, middle)

if ~exist('cmap', 'var') || isempty(cmap)
	cmap = 'spring';
end

if ~exist('middle', 'var') || isempty(middle)
	middle = 1;
end

map = eval([cmap '(3)']);
map(2,:) = middle.*[1 1 1];
tempG = interp1([1 20], map(1:2, 2), 1:20);
tempB = interp1([1 20], map(1:2, 3), 1:20);
tempG2 = interp1([1 20], map(2:3, 2), 1:20);
tempB2 = interp1([1 20], map(2:3, 3), 1:20);
tempR = interp1([1 20], map(1:2, 1), 1:20);
tempR2 = interp1([1 20], map(2:3, 1), 1:20);
map = [tempR' tempG' tempB'; tempR2' tempG2' tempB2'];
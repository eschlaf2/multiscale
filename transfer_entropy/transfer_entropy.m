% ?Lee, Jongsun, Hyun-Soo Choi, Yongkweon Jeon, Yongsik Kwon, Donghun Lee,
% and Sungroh Yoon.  2018.  ?Detecting System Anomalies in Multivariate
% Time Series with Information Transfer and Random Walk?

% ?Lizier, Joseph T.  2014.  ?JIDT: An Information-Theoretic Toolkit for
% Studying the Dynamics of Complex Systems?

% Change location of jar to match yours:
addpath(genpath('/Users/emilyschlafly/Documents/GITHUB/jidt/'));
javaaddpath('/Users/emilyschlafly/Documents/GITHUB/jidt/infodynamics.jar')

%%
TSTART = 0;  % s
TEND = Inf;  % s
FREQ = mea.SamplingRate;  % Hz
WINDOW = 1;  % s
SHIFT = .1;  % s
DATA = 'mua';  % which field to use from mea
C = 0.05;  % restart probability

%%
% Convert parameter units to samples
window = WINDOW * FREQ;
shift = SHIFT * FREQ;
t = TSTART * FREQ + 1 : shift : min(TEND * FREQ, size(mea.(DATA), 1) - window);
numCh = size(mea.(DATA), 2);
q = numel(t);

teCalc=javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorMultiVariateKraskov');
teCalc.initialise(1,sourceDim,destDim); % Use history length 1 (Schreiber k=1)
teCalc.setProperty('k', '4'); % Use Kraskov parameter K=4 for 4 nearest points

res = zeros(numel(t), numCh, numCh);
R = res;
similarity = zeros(numel(t) - 1, 1);

i = 1;
for tt = t
	obsn = tt : tt + window - 1;
	for chS = 1:numCh
		for chD = 1:numCh
% 			if chD == chS
% 				continue
% 			end
			sourceCh = chS;
			destCh = chD;

			sourceDim = numel(sourceCh);
			destDim = numel(destCh);			

			% Generate some random normalised data.
			% % numObservations = 10000;
		% 	covariance=0.4;

			window = length(obsn);  % EDS 1/16

			% % sourceMVArray = randn(numObservations, sourceDim);
			sourceMVArray = mea.(DATA)(obsn, sourceCh);  % EDS 1/16

			% Set first two columns of dest to copy source values
			% % destMVArray  = [zeros(1,sourceDim); covariance*(sourceMVArray(1:numObservations-1,:)) + (1-covariance)*randn(numObservations-1, sourceDim)];
			% Set a third colum to be randomised
			% % destMVArray(:,3) = randn(numObservations, 1);
			% % sourceMVArray2= randn(numObservations, sourceDim); % Uncorrelated source

			destMVArray = mea.Data(obsn, destCh);

			% Create a TE calculator and run it:
			teCalc.setObservations(octaveToJavaDoubleMatrix(sourceMVArray), octaveToJavaDoubleMatrix(destMVArray));
			% Perform calculation with correlated source:
			res(i, chS, chD) = teCalc.computeAverageLocalOfObservations();
			% Note that the calculation is a random variable (because the generated
			%  data is a set of random variables) - the result will be of the order
			%  of what we expect, but not exactly equal to it; in fact, there will
			%  be some variance around it. It will probably be biased down here
			%  due to small correlations between the supposedly uncorrelated variables.
		% 	fprintf('TE result %.4f nats; expected to be close to %.4f nats for the two correlated Gaussians\n', ...
		% 	result, 2*log(1/(1-covariance^2)));

			% Perform calculation with uncorrelated source:
		% 	teCalc.initialise(1,sourceDim,destDim); % Initialise leaving the parameters the same
		% 	teCalc.setObservations(octaveToJavaDoubleMatrix(sourceMVArray2), octaveToJavaDoubleMatrix(destMVArray));
		% 	result2 = teCalc.computeAverageLocalOfObservations();
		% 	fprintf('TE result %.4f nats; expected to be close to 0 nats for these uncorrelated Gaussians\n', result2);
		end
		res(i, chS, :) = res(i, chS, :) / (log(q) + sum(res(i, chS, :), 3));
	end
	temp = C * inv(eye(numCh) - (1 - C) * squeeze(res(i, :, :)));
	R(i, :, :) = temp - diag(diag(temp));
	if i > 1
		similarity(i - 1) = 1 / (1 + norm(squeeze(R(i, :, :) - R(i - 1, :, :)), 'fro'));
	end
	i = i + 1;
end
clear teCalc
% plot((t / FREQ) - mea.Padding(1), res);
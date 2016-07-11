function dragdata = readdata(species)

% Read in the data
fname = dir(['data/Training Set/*_',species,'_*']);
drag = load(['data/Training Set/',fname.name]);

% Get the output;
ysim = drag(:,end)';
ysimmean = mean(ysim);
ysimStd = ysim - ysimmean;
ysimsd = std(ysimStd);
ysimStd = ysimStd / ysimsd;


x = drag(:,1:(end-1));
xmin = min(x);
xrange = range(x);
x = bsxfun(@minus, x, xmin);
x = bsxfun(@rdivide, x, xrange);

simData.x = x;
simData.yStd = ysimStd;
simData.Ksim = 1;
simData.orig.ysim = ysim;
simData.orig.xmin = xmin;
simData.orig.xrange = xrange;
simData.orig.ymean = ysimmean;
simData.orig.ysd = ysimsd;

dragdata = struct('simData',simData,'obsData',[]);
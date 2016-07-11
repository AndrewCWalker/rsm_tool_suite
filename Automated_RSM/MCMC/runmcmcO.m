species = 'O';

ddat = readdata(species);
params = setupModel(ddat.obsData, ddat.simData);

% Fix lamWs prior
params.priors.lamWs.params = [ones(params.model.pu,1), zeros(params.model.pu,1)];
params.priors.lamWs.bUpper = Inf;
params.priors.lamWs.bLower = 0;

% Fix lamWOs
params.mcmc.svars = params.mcmc.svars(1:3);
params.mcmc.svarSize = params.mcmc.svarSize(1:3);
params.mcmc.wvars = params.mcmc.wvars(1:3);
params.model.lamWOs = 1e9;

% Step-size estimation, burn-in, estimation
nburn =100; nlev = 25; nmcmc = 10000;
params = stepsize(params,nburn,nlev);
pout = gpmmcmc(params,nmcmc,'step',1);
save(['pout/pout',species,'.mat'], 'pout');

%load(['pout',species,'.mat'], 'pout');
nburn =100; nlev = 25; nmcmc = 10000;
pvec = floor(linspace(nburn*nlev+1, nburn*nlev+nmcmc, 500));
pu = pout.model.pu;
m = pout.model.m;
p = pout.model.p;
xmin = pout.simData.orig.xmin;
xrange = pout.simData.orig.xrange;

% Make some boxplots of the correlations.
rhoU = exp(-[pout.pvals((nburn*nlev+1):(nburn*nlev+nmcmc)).betaU]/4)';
for i=1:pu
    subplot(pu,1,i);
    ix = (i-1)*p + (1:p);
    boxplot(rhoU(:,ix));
    ylim([0,1]);
end
print('-depsc', ['rho/rho',species,'.eps']);
close;

% % Test predictions.
% Take the mean of pout
pout = meanpout(pout, pvec);
pvec = 1;

fname = dir(['data/Test Set/*_',species,'_*']);
testdata = load(['data/Test Set/',fname.name]);
testy = testdata(:,end);
testx = testdata(:,1:(end-1));
testx = bsxfun(@minus, testx, xmin);
testx = bsxfun(@rdivide, testx, xrange);
mtest = length(testy);

testpred = zeros(mtest,1);
for i=1:mtest
    disp(i);
    wpred = gPredict(testx(i,:), pout.pvals(pvec), pout.model, pout.data, 'returnMuSigma',1);
    testpred(i) = wpred.Myhat*pout.simData.orig.ysd + pout.simData.orig.ymean;
end
save(['testpred/testpred',species,'.mat'], 'testpred');

plot(testy, testpred, 'o');
xl = xlim;
yl = ylim;
lim = [min(xl(1), yl(2)), max(xl(2), yl(2))];
hold on
plot(lim, lim, '--r');
title('Test Set Predictions', 'FontSize',16);
xlabel('Data', 'FontSize', 14)
ylabel('Emulator', 'FontSize', 14)
print('-depsc', ['testpred/testpred',species,'.eps']);
close;

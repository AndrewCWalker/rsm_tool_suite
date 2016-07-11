function poutmean = meanpout(pout, pvec)

poutmean = pout;
pvals.betaU = mean([pout.pvals(pvec).betaU],2);
pvals.lamUz = mean([pout.pvals(pvec).lamUz],2);
pvals.lamWs = mean([pout.pvals(pvec).lamWs],2);
pvals.lamWOs = mean([pout.pvals(pvec).lamWOs],2);


poutmean.pvals = pvals;

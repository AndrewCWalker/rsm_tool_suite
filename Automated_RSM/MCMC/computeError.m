speciesList = {'H', 'He', 'N', 'N2', 'O', 'O2'};

rmse = zeros(length(speciesList),1);
rmspe = zeros(length(speciesList),1);
for i=1:length(speciesList)
    species = speciesList{i};
    fname = dir(['data/Test Set/*_',species,'_*']);
    testdata = load(['data/Test Set/',fname.name]);
    testpred = load(['testpred/testpred',species,'.mat']);
    res = testpred.testpred - testdata(:,end);
    rmse(i) = sqrt(mean(res.^2));
    rmspe(i) = 100*sqrt(mean((res./testdata(:,end)).^2));
end
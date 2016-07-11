species = {'H', 'He', 'N', 'N2', 'O', 'O2'};

p = 7;
pvec = 2501:12500;

% Gather some stuff together
M = 100;
x = zeros(M, p, length(species));
xmin = zeros(length(species), p);
xrange = zeros(length(species), p);
corrlength = zeros(length(species),p);
meansd = zeros(2,length(species));
precision = zeros(2, length(species));
stddata = zeros(length(species), M);
mm = zeros(length(species),1);

% The parameter order needs to be rearranged.
porder = 1:7;

for i=1:length(species)
    % Load the GPMSA object.
    load(['pout/pout',species{i},'.mat'], 'pout');
    
    % Size
    m = pout.model.m;
    mm(i) = m;
    
    % Get the min and max of the design.
    xmin(i,:) = pout.simData.orig.xmin(porder);
    xrange(i,:) = pout.simData.orig.xrange(porder);
    
    % Get the design
    x(1:m,:,i) = pout.simData.x(:,porder);
    
    % Get the correlation lengths
    tmp = mean([pout.pvals(pvec).betaU],2)';
    corrlength(i,:) = tmp(porder);
    
    % Precisions
    precision(1,i) = mean([pout.pvals(pvec).lamUz]);
    precision(2,i) = mean(1./(1./[pout.pvals(pvec).lamWs] + 1./[pout.pvals(pvec).lamWOs]));
    
    % Standardized training data
    stddata(i,1:m) = pout.data.w';
    
    % Means and sds
    meansd(1,i) = pout.simData.orig.ymean;
    meansd(2,i) = pout.simData.orig.ysd;
end

%%%%%%% Write the RSM file %%%%%%%

%Write header
txt = ['#RSM Parameters'];
dlmwrite('RSM.dat', txt, 'delimiter', '')
txt = [' '];
dlmwrite('RSM.dat', txt, '-append')

% beta
txt = ['#beta[',num2str(length(species)),'][',num2str(p),']'];
dlmwrite('RSM.dat',txt, '-append', 'delimiter', '');
tmp = [corrlength];
dlmwrite('RSM.dat', tmp, '-append', 'delimiter', ' ', 'precision', '%1.8e');

% xmin
txt = [' '];
dlmwrite('RSM.dat', txt, '-append')
txt = ['#xmin[',num2str(length(species)),'][',num2str(p),']'];
dlmwrite('RSM.dat',txt,'-append','delimiter','');
tmp = [xmin];
dlmwrite('RSM.dat', tmp,'-append', 'delimiter', ' ', 'precision', '%1.8e');

% xrange
txt = [' '];
dlmwrite('RSM.dat', txt, '-append')
txt = ['#xrange[',num2str(length(species)),'][',num2str(p),']'];
dlmwrite('RSM.dat',txt,'-append','delimiter','');
tmp = [xrange];
dlmwrite('RSM.dat', tmp,'-append', 'delimiter', ' ', 'precision', '%1.8e');

% Means
txt = [' '];
dlmwrite('RSM.dat', txt, '-append')
txt = ['#mean[',num2str(length(species)),']'];
dlmwrite('RSM.dat',txt,'-append','delimiter','');
tmp = [meansd(1,:)];
dlmwrite('RSM.dat', tmp,'-append', 'delimiter', ' ', 'precision', '%1.8e');

% SDs
txt = [' '];
dlmwrite('RSM.dat', txt, '-append')
txt = ['#sd[',num2str(length(species)),']'];
dlmwrite('RSM.dat', txt, '-append', 'delimiter', '');
tmp = [meansd(2,:)];
dlmwrite('RSM.dat', tmp, '-append', 'delimiter', ' ', 'precision', '%1.8e');

% lamz
txt = [' '];
dlmwrite('RSM.dat', txt, '-append')
txt = ['#lamz[',num2str(length(species)),']'];
dlmwrite('RSM.dat',txt,'-append','delimiter','');
tmp = [precision(1,:)];
dlmwrite('RSM.dat', tmp,'-append', 'delimiter', ' ', 'precision', '%1.8e');

% lamws
txt = [' '];
dlmwrite('RSM.dat', txt, '-append')
txt = ['#lamws[',num2str(length(species)),']'];
dlmwrite('RSM.dat', txt, '-append', 'delimiter', '');
tmp = [precision(2,:)];
dlmwrite('RSM.dat', tmp, '-append', 'delimiter', ' ', 'precision', '%1.8e');

% x
txt = [' '];
dlmwrite('RSM.dat', txt, '-append')
txt = ['#x[',num2str(length(species)),'][',num2str(M),'][', num2str(p),']'];
dlmwrite('RSM.dat', txt, '-append', 'delimiter', '');
for i=1:length(species)
    if i>1
        txt = [' '];
        dlmwrite('RSM.dat', txt', '-append')
    end
    txt = ['#',species{i}];
    dlmwrite('RSM.dat',txt,'-append', 'delimiter', '');
    tmp = [x(:,:,i)];
    dlmwrite('RSM.dat', tmp, '-append', 'delimiter', ' ', 'precision', '%1.8e');
end

% w
txt = [' '];
dlmwrite('RSM.dat', txt, '-append')
txt = ['#wt[',num2str(M),'][',num2str(length(species)),']'];
dlmwrite('RSM.dat', txt, '-append', 'delimiter', '');
tmp = [transpose(stddata)];
dlmwrite('RSM.dat', tmp,'-append', 'delimiter', ' ', 'precision', '%1.8e');




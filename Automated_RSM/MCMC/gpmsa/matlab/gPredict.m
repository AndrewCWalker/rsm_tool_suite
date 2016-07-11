%function pred=gPredict(xpred,pvals,model,data,varargs)
% Predict using a gpmsa constructed model. 
% result is a 3-dimensional prediction matrix:
%       #pvals by model-dims by length-xpred 
%   model-dims is the simulation basis size for a w-prediction (a model
%   with no observation data) or the v (discrepancy basis loadings) and u 
%   (simulation basis loadings) for uv-prediction
% argument pair keywords are
%   mode - 'wpred' or 'uvpred'
%   theta - the calibration variable value(s) to use, either one value or
%           one for each xpred. Default is the calibrated param in pvals
%   addResidVar - default 0, whether to add the residual variability
%   returnRealization - default 1, whether to return realizations of the
%           process specified. values will be in a .w field for wpred, or
%           .v and .u fields for uvpred.
%   returnMuSigma - default 0, whether to return the mean and covariance of
%           the process specified. results will be in a matrix .Myhat field
%           and a cell array .Syhat field. 
%   returnCC - default 0, return the cross-covariance of the data and
%           predictors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: James R. Gattiker, Los Alamos National Laboratory
%
% This file was distributed as part of the GPM/SA software package
% Los Alamos Computer Code release LA-CC-06-079, C-06,114
%
% Copyright 2008.  Los Alamos National Security, LLC. This material 
% was produced under U.S. Government contract DE-AC52-06NA25396 for 
% Los Alamos National Laboratory (LANL), which is operated by Los Alamos 
% National Security, LLC for the U.S. Department of Energy. The U.S. 
% Government has rights to use, reproduce, and distribute this software.  
% NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY 
% WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF 
% THIS SOFTWARE.  If software is modified to produce derivative works, 
% such modified software should be clearly marked, so as not to confuse 
% it with the version available from LANL.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pred=gPredict(xpred,pvals,model,data,varargin)

  if model.n==0
    mode='wpred';
  else
    mode='uvpred';
  end
  theta=[];
  addResidVar=0;
  returnRealization=1;
  returnMuSigma=0;
  returnCC=0;
  parseAssignVarargs({'mode','theta','addResidVar', ...
                      'returnRealization','returnMuSigma', ...
                      'returnCC'}); 

  switch mode
    case 'uvpred'
        pred=uvPred(xpred,pvals,model,data,theta,addResidVar, ...
                    returnRealization,returnMuSigma,returnCC);
    case 'wpred'
        pred=wPred(xpred,pvals,model,data,theta,addResidVar, ...
                    returnRealization,returnMuSigma,returnCC);
    otherwise
      error('invalid mode in gPredict');
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pred=wPred(xpred,pvals,model,data,thetapred,addResidVar,retRlz,retMS,retC)

  n=model.n; m=model.m;p=model.p;q=model.q;pu=model.pu;

  npred=size(xpred,1);

  diags1=diagInds(m*pu);
  diags2=diagInds(npred*pu);

  nreal=length(pvals);
  tpred=zeros([nreal,npred*pu]);

  for ii=1:length(pvals)
    if n>0
      theta=pvals(ii).theta'; 
    end
    betaU=reshape(pvals(ii).betaU,p+q,pu);
    lamUz=pvals(ii).lamUz;
    lamWs=pvals(ii).lamWs; lamWOs=pvals(ii).lamWOs;

    if n>0
      if isempty(thetapred)
        xpredt=[xpred repmat(theta,npred,1)];
      else
        xpredt=[xpred thetapred];
      end
    else
      xpredt=xpred;
    end

    xpredDist=genDist(xpredt);  
    zxpredDist=genDist2(data.zt,xpredt);

      SigW=zeros(m*pu);
        for jj=1:pu
          bStart=(jj-1)*m+1; bEnd=bStart+m-1; 
          SigW(bStart:bEnd,bStart:bEnd)=...
              gCovMat(model.ztDist,betaU(:,jj),lamUz(jj));
        end
        SigW(diags1)=SigW(diags1)+ ...
            kron(1./(model.LamSim*lamWOs)',ones(1,m)) + ...
            kron(1./(lamWs)',ones(1,m)) ;

      SigWp=zeros(npred*pu);
        for jj=1:pu
          bStart=(jj-1)*npred+1; bEnd=bStart+npred-1;
          SigWp(bStart:bEnd,bStart:bEnd)= ...
              gCovMat(xpredDist,betaU(:,jj),lamUz(jj));
        end
        SigWp(diags2)=SigWp(diags2)+ ...
             kron(1./(lamWs)',ones(1,npred)) ;
        if addResidVar
          SigWp(diags2)=SigWp(diags2)+ ...
             kron(1./(model.LamSim*lamWOs)',ones(1,npred));
        end
      SigWWp=zeros(m*pu,npred*pu);
        for jj=1:pu
          bStartI=(jj-1)*m+1; bEndI=bStartI+m-1;
          bStartJ=(jj-1)*npred+1; bEndJ=bStartJ+npred-1;
          SigWWp(bStartI:bEndI,bStartJ:bEndJ)=...
              gCovMat(zxpredDist,betaU(:,jj),lamUz(jj));
        end

    SigData=SigW;
    SigPred=SigWp;
    SigCross=SigWWp;
    % Get the stats for the prediction stuff. 
      %W=(SigCross')/SigData;
      W=linsolve(SigData,SigCross,struct('SYM',true,'POSDEF',true))';
      Myhat=W*(data.w(:));
      Syhat=SigPred-W*SigCross;
      
    if retRlz
      % And do a realization
      tpred(ii,:)=rmultnormsvd(1,Myhat,Syhat')';
    end
    if retMS
      % add the distribution params
      pred.Myhat(ii,:)=Myhat;
      pred.Syhat{ii}=Syhat;
    end
    if retC
      pred.CC{ii}=SigCross;
    end
  end
  
  if retRlz
    % Reshape the pred matrix to 3D:
    %  first dim  - (number of realizations [pvals])
    %  second dim - (number of principal components)
    %  third dim  - (number of points [x,theta]s) 
    pred.w=zeros(length(pvals),pu,npred);
    for ii=1:pu
      pred.w(:,ii,:)=tpred(:,(ii-1)*npred+1:ii*npred);
    end
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pred=uvPred(xpred,pvals,model,data,thetapred,addResidVar,retRlz,retMS,retC)
  n=model.n;m=model.m;p=model.p;q=model.q;pu=model.pu;pv=model.pv;
  lamVzGnum=model.lamVzGnum; lamVzGroup=model.lamVzGroup;

  npred=size(xpred,1);

  diags0=diagInds(n*pu);
  diags1=diagInds(m*pu);
  diags2=diagInds(npred*pu);

  x0Dist=genDist(data.x);
  xpred0Dist=genDist(xpred);
  xxpred0Dist=genDist2(data.x,xpred);

  nreal=length(pvals);
  tpred=zeros([nreal,npred*(pv+pu)]);

  for ii=1:length(pvals)
    theta=pvals(ii).theta';
    betaV=reshape(pvals(ii).betaV,p,lamVzGnum); 
    betaU=reshape(pvals(ii).betaU,p+q,pu);
    lamVz=pvals(ii).lamVz; lamUz=pvals(ii).lamUz; lamWOs=pvals(ii).lamWOs;
    lamWs=pvals(ii).lamWs; lamOs=pvals(ii).lamOs;

    if isempty(thetapred)
      xpredt=[xpred repmat(theta,npred,1)];
    else
      xpredt=[xpred thetapred];
    end

    xDist=genDist([data.x repmat(theta,n,1)]);
    ztDist=genDist(data.zt);
    xzDist=genDist2([data.x repmat(theta,n,1)],data.zt);
    xpredDist=genDist(xpredt);  
    xxpredDist=genDist2([data.x repmat(theta,n,1)],xpredt);
    zxpredDist=genDist2(data.zt,xpredt);

    % Generate the part of the matrix related to the data
    % Four parts to compute: Sig_v, Sig_u, Sig_w, and the Sig_uw crossterm
      SigV=zeros(n*pv);
        for jj=1:lamVzGnum;
          vCov(jj).mat=gCovMat(x0Dist, betaV(:,jj), lamVz(jj));
        end
        for jj=1:pv
          bStart=(jj-1)*n+1; bEnd=bStart+n-1;
          SigV(bStart:bEnd,bStart:bEnd)=vCov(lamVzGroup(jj)).mat;
        end
      SigU=zeros(n*pu);
        for jj=1:pu
          bStart=(jj-1)*n+1; bEnd=bStart+n-1;
          SigU(bStart:bEnd,bStart:bEnd)= ...
              gCovMat(xDist,betaU(:,jj),lamUz(jj));
        end
        SigU(diags0)=SigU(diags0)+...
             kron(1./(lamWs)',ones(1,n)) ;
      SigW=zeros(m*pu);
        for jj=1:pu
          bStart=(jj-1)*m+1; bEnd=bStart+m-1; 
          SigW(bStart:bEnd,bStart:bEnd)=...
              gCovMat(ztDist,betaU(:,jj),lamUz(jj));
        end
        SigW(diags1)=SigW(diags1)+ ...
            kron(1./(model.LamSim*lamWOs)',ones(1,m)) + ...
            kron(1./(lamWs)',ones(1,m)) ;
      SigUW=zeros(n*pu,m*pu);
        for jj=1:pu
          bStartI=(jj-1)*n+1; bEndI=bStartI+n-1;
          bStartJ=(jj-1)*m+1; bEndJ=bStartJ+m-1;
          SigUW(bStartI:bEndI,bStartJ:bEndJ)=...
              gCovMat(xzDist,betaU(:,jj),lamUz(jj));
        end
        if model.scOut
          SigData=[ SigU+SigV    SigUW; ...
                    SigUW'       SigW ];
          SigData(1:n*pu,1:n*pu) = ...
            SigData(1:n*pu,1:n*pu) + model.SigObs*1/lamOs;
        else
          SigData=[SigV                 zeros(n*pv,(n+m)*pu);  ...
                   zeros((n+m)*pu,n*pv) [ SigU    SigUW; ...
                                          SigUW'  SigW  ] ];
          SigData(1:n*(pv+pu),1:n*(pv+pu)) = ...
            SigData(1:n*(pv+pu),1:n*(pv+pu)) + model.SigObs*1/lamOs;
        end

    % Generate the part of the matrix related to the predictors
    % Parts to compute: Sig_vpred, Sig_upred
      SigVp=zeros(npred*pv);
        for jj=1:lamVzGnum;
          vpCov(jj).mat=gCovMat(xpred0Dist, betaV(:,jj), lamVz(jj));
        end
        for jj=1:pv
          bStart=(jj-1)*npred+1; bEnd=bStart+npred-1;
          SigVp(bStart:bEnd,bStart:bEnd)=vpCov(lamVzGroup(jj)).mat;
        end
        %SigVp(diagInds(npred*pv))=SigVp(diagInds(npred*pv))+1;
      SigUp=zeros(npred*pu);
        for jj=1:pu
          bStart=(jj-1)*npred+1; bEnd=bStart+npred-1;
          SigUp(bStart:bEnd,bStart:bEnd)= ...
              gCovMat(xpredDist,betaU(:,jj),lamUz(jj));
        end
        SigUp(diags2)=SigUp(diags2)+...
             kron(1./(lamWs)',ones(1,npred)) ;
        if addResidVar
          SigUp(diags2)=SigUp(diags2)+ ...
             kron(1./(model.LamSim*lamWOs)',ones(1,npred)) ;
        end
           
      SigPred=[SigVp                     zeros(npred*pv,npred*pu);  ...
               zeros(npred*pu,npred*pv)  SigUp  ];

    % Now the cross-terms. 
      SigVVx=zeros(n*pv,npred*pv);    
        for jj=1:lamVzGnum;
          vvCov(jj).mat=gCovMat(xxpred0Dist, betaV(:,jj), lamVz(jj));
        end
        for jj=1:pv
          bStartI=(jj-1)*n+1; bEndI=bStartI+n-1;
          bStartJ=(jj-1)*npred+1; bEndJ=bStartJ+npred-1;
          SigVVx(bStartI:bEndI,bStartJ:bEndJ)=vvCov(lamVzGroup(jj)).mat; 
        end
      SigUUx=zeros(n*pu,npred*pu);
        for jj=1:pu
          bStartI=(jj-1)*n+1; bEndI=bStartI+n-1;
          bStartJ=(jj-1)*npred+1; bEndJ=bStartJ+npred-1;
          SigUUx(bStartI:bEndI,bStartJ:bEndJ)=...
              gCovMat(xxpredDist,betaU(:,jj),lamUz(jj));
        end
      SigWUx=zeros(m*pu,npred*pu);
        for jj=1:pu
          bStartI=(jj-1)*m+1; bEndI=bStartI+m-1;
          bStartJ=(jj-1)*npred+1; bEndJ=bStartJ+npred-1;
          SigWUx(bStartI:bEndI,bStartJ:bEndJ)=...
              gCovMat(zxpredDist,betaU(:,jj),lamUz(jj));
        end
      if model.scOut
        SigCross=[SigVVx                 SigUUx; ...
                  zeros(m*pu,npred*pv)   SigWUx];
      else
        SigCross=[SigVVx                 zeros(n*pv,npred*pu); ...
                  zeros(n*pu,npred*pv)   SigUUx; ...
                  zeros(m*pu,npred*pv)   SigWUx];
      end

    if 0
      figure(3)
      subplot(2,2,1); imagesc(gScale(SigData,'sqrt'))
      subplot(2,2,2); imagesc(gScale(SigCross,'sqrt'))
      subplot(2,2,3); imagesc(gScale(SigCross','sqrt'))
      subplot(2,2,4); imagesc(gScale(SigPred,'sqrt'))
      keyboard
    end

    % Get the stats for the prediction stuff. 
      %W=(SigCross')/SigData;
      W=linsolve(SigData,SigCross,struct('SYM',true,'POSDEF',true))';
      if model.scOut, Myhat=W*model.uw; else Myhat=W*model.vuw; end
      Syhat=SigPred-W*SigCross;
      
    if retRlz
      % And do a realization
      tpred(ii,:)=rmultnormsvd(1,Myhat,Syhat')';
    end
    if retMS
      % log the distribution params
      pred.Myhat(ii,:)=Myhat;
      pred.Syhat{ii}=Syhat;
    end
    if retC
      pred.CC{ii}=SigCross;
    end
  end

  if retRlz
    % Reshape the pred matrix to 3D, for each component:
    %  first dim  - (number of realizations [pvals])
    %  second dim - (number of principal components)
    %  third dim  - (number of points [x,theta]s) 
    pred.v=zeros(length(pvals),pv,npred);
    pred.u=zeros(length(pvals),pu,npred);
    for ii=1:pv
      pred.v(:,ii,:)=tpred(:,(ii-1)*npred+1:ii*npred);
    end
    for ii=1:pu
      pred.u(:,ii,:)=tpred(:,pv*npred+((ii-1)*npred+1:ii*npred) );
    end
  end

end


%
% Helper function rmultnormSVD computes multivariate normal realizations
function rnorm = rmultnormsvd(n,mu,cov)
  [U S] = svd(cov);
  rnorm = repmat(mu,1,n) + U*sqrt(S) * randn(size(mu,1),n);
end 

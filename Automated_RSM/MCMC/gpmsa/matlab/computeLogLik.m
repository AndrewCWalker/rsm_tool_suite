% function model = computeLogLik(model,data,C)
% 
% Builds the log likelihood of the data given the model parameters.
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = computeLogLik(model,data,C)

 n=model.n;   m=model.m;
pu=model.pu; pv=model.pv;
 p=model.p;   q=model.q;
lamVzGnum=model.lamVzGnum; lamVzGroup=model.lamVzGroup;

% validate and process the changed field
do.theta=0;do.betaV=0;do.lamVz=0;do.betaU=0;do.lamUz=0;
do.lamWs=0;do.lamWOs=0;

if strcmp(C.var,'theta') || strcmp(C.var,'all') % update distances
       model.xDist=genDist([data.x repmat(model.theta,n,1)]);
       model.xzDist=genDist2([data.x repmat(model.theta,n,1)],data.zt);
end

switch C.var
    case 'all'
       do.theta=1; do.betaV=1; do.lamVz=1; do.betaU=1; 
       do.lamUz=1; do.lamWs=1; do.lamWOs=1; do.lamOs=1; 
       model.SigWl=zeros(pu,1);       
    case 'theta'; do.theta=1;
    case 'betaV'; do.betaV=1;
    case 'lamVz'; do.lamVz=1;
    case 'betaU'; do.betaU=1;
    case 'lamUz'; do.lamUz=1;
    case 'lamWs'; do.lamWs=1;
    case 'lamWOs'; do.lamWOs=1;
    case 'lamOs'; do.lamOs=1;
    otherwise
        %error('Invalid Subtype in computeLogLik');
end

betaV=model.betaV;   lamVz=model.lamVz;
betaU=model.betaU;   lamUz=model.lamUz;
lamWs=model.lamWs;   lamWOs=model.lamWOs; lamOs=model.lamOs;

% Four parts to compute: Sig_v, Sig_u, Sig_w, and the Sig_uw crossterm

  if (do.theta || do.betaV || do.lamVz)
	SigV=[];
    for jj=1:lamVzGnum
      SigV(jj).mat=gCovMat(model.x0Dist, betaV(:,jj), lamVz(jj));
    end
    model.SigV=SigV;
  else
    SigV=model.SigV;
  end

  if (do.theta || do.betaU || do.lamUz || do.lamWs)
    SigU(pu).mat=[];
    diags1=diagInds(n);
    for jj=1:pu
        SigU(jj).mat=gCovMat(model.xDist,betaU(:,jj),lamUz(jj));
        SigU(jj).mat(diags1)=SigU(jj).mat(diags1)+1/lamWs(jj);
    end
    model.SigU=SigU;
  else
    SigU=model.SigU;
  end

  if (do.betaU || do.lamUz || do.lamWs || do.lamWOs)
    diags2=diagInds(m);
    switch C.var
      case 'all'
        jinds=1:pu;
      case 'betaU'
        jinds=ceil( C.index/(p+q) );
      case {'lamUz','lamWs'}
        jinds=C.index;
      case 'lamWOs'
        jinds=1:pu;
    end
    for jj=jinds
      if isempty(data.ztSep)
        cg=gCovMat(model.ztDist,betaU(:,jj),lamUz(jj));
        cg(diags2)=cg(diags2)+1/(model.LamSim(jj)*lamWOs) + 1/lamWs(jj);

        % calculate the SigW likelihood for each block
          model.SigWl(jj)=doLogLik(cg,model.w((jj-1)*m+1:jj*m));
        % calculate the SigW inverse for each block 
          model.SigWi(jj).mat=inv(cg);
      else
        % there is a separable design, so compute these as kron'ed blocks
        segVarStart=1;
        for ii=1:length(data.ztSep)
           segVars=segVarStart:(segVarStart + model.ztSepDist{ii}.p-1);
           segVarStart=segVarStart+ model.ztSepDist{ii}.p;
           cg{ii}=gCovMat(model.ztSepDist{ii},betaU(segVars,jj),lamUz(jj));
        end
        cgNugget=1/(model.LamSim(jj)*lamWOs) + 1/lamWs(jj);
        [model.SigWl(jj) model.SigWi(jj).mat]= ...
            doLogLikSep(cg,cgNugget,model.w((jj-1)*m+1:jj*m));
      end
    end
  end

  if (do.theta || do.betaU || do.lamUz)
    SigUW(pu).mat=[];
    for jj=1:pu
        SigUW(jj).mat=gCovMat(model.xzDist,betaU(:,jj),lamUz(jj));
    end
    model.SigUW=SigUW;
  else
    SigUW=model.SigUW;
  end
  
  % The computation is decomposed into the likelihood of W, 
  %  and the likelihood of VU|W. 

  % Compute the likelihood of the W part (already done the blocks)
    LogLikW=sum(model.SigWl);

  % Compute the likelihood of the VU|W
  % This requires using the gaussian model estimation stuff.
    % It is complicated because of shortcuts allowed by lack of correlation
    % between W and V

    % do these ops, on the block diagonal blocks: 
    %    W=SigUW*model.SigWi;
    %    SigUgW=SigU-W*SigUW';
      W(pu).mat=[];
      SigUgW(pu).mat=[];
      for ii=1:pu
        W(ii).mat=SigUW(ii).mat*model.SigWi(ii).mat;
        SigUgW(ii).mat=SigU(ii).mat-W(ii).mat*SigUW(ii).mat';
      end

    %for scalar output: SigVUgW=[SigV+SigUgW] ...
    %                          + model.SigObs/lamOs;
    %otherwise:         SigVUgW=[SigV             zeros(n*pv,n*pu); ...
    %                            zeros(n*pu,n*pv) SigUgW         ] ...
    %                          + model.SigObs/lamOs;
      SigVUgW=model.SigObs/lamOs;
        for ii=1:pv
          blkRange=(ii-1)*n+1:ii*n;
          SigVUgW(blkRange,blkRange)=SigVUgW(blkRange,blkRange)+ ...
                                                  SigV(lamVzGroup(ii)).mat;
        end
      if model.scOut
        for ii=1:pu
          blkRange=(ii-1)*n+1:ii*n;
          SigVUgW(blkRange,blkRange)=SigVUgW(blkRange,blkRange)+SigUgW(ii).mat;
        end
      else
        for ii=1:pu
          blkRange=n*pv+[(ii-1)*n+1:ii*n];
          SigVUgW(blkRange,blkRange)=SigVUgW(blkRange,blkRange)+SigUgW(ii).mat;
        end
      end

    % do this op: MuVUgW =W*model.w;
      MuVUgW=zeros(n*pu,1);
      for ii=1:pu
        blkRange1=(ii-1)*n+1:ii*n;
        blkRange2=(ii-1)*m+1:ii*m;
        MuVUgW(blkRange1)=W(ii).mat*model.w(blkRange2);
      end
    
    % for scalar output:  MuDiff=   [u] - [MuVUgW]
    % otherwise:          MuDiff= [v;u] - [0;MuVUgW] 
      if model.scOut
         MuDiff=model.u; MuDiff=MuDiff-MuVUgW;
      else
         MuDiff=model.vu; MuDiff(pv*n+1:end)=MuDiff(pv*n+1:end)-MuVUgW;
      end

    % Now we can get the LL of VU|W
      LogLikVUgW=doLogLik(SigVUgW,MuDiff);
  
  % Final Answer, LL(VU) = LL(VU|W) + LL(W)
    model.logLik=LogLikVUgW+LogLikW;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L chCov]=doLogLik(Sigma,data)
  chCov=chol(Sigma);
  logDet=sum(log(diag(chCov))); % actually, the log sqrt(det)
  p1=(chCov')\data;
  L=-logDet-0.5*(p1'*p1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L Sinv]=doLogLikSep(Sigma,nugget,data)

  % get the eigenDecomp for the blocks
    for ii=1:length(Sigma)
      [V{ii} D{ii}]=eig(Sigma{ii});
    end
  
  % compute the determinant from the eigenvalues and the nugget
    dkron=diag(D{end});
    for ii=(length(Sigma)-1):-1:1
      dkron=kron(diag(D{ii}),dkron);
    end
    logDet=log(prod(dkron+nugget));

  % compute the composite inverse, including the nugget
    vkron=V{end};
    for ii=(length(Sigma)-1):-1:1
      vkron=kron(V{ii},vkron);
    end
    Sinv=vkron * diag(1./(dkron+nugget)) * vkron';

  % compute the log likelihood
    L=-0.5*logDet-0.5*data'*Sinv*data;

end

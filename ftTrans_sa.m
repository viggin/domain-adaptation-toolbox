function [ftAllNew,transMdl] = ftTrans_sa(ft,maSrc,target,maLabeled,param)
%Subspace Alignment (SA)
% 
%	A feaure-level transfer learning (domain adaptation) algorithm which
% matches two subspaces. Application scope:
%	+ two discrete domains (source and target)
%	+ unlabeled source domain
%	+ unlabeled target domain
%
% ftAll:	All samples in all domains. n-by-m matrix, n is the number of 
%	samples, m is the dimension of features.
% maSrc:	n-by-1 logical vector, maSrc(i)=true if i is from the source
%	domain, 0 if target domain.
% target and maLabeled: useless
% 
% param: Struct of hyper-parameters, please see the first cell of this
%	program ("default parameters") for details. You can set parameter p to 
%	x by setting param.p = x. For parameters that are not set, default 
%	values will be used.
% 
% ftAllNew:	All samples in the learned subspace.
% transMdl:	A struct containing the model, transMdl.W is the projection
%	matrix.
% 
% The ftProc_pca_tr function is a part of the PRTools toolbox.
% 
% ref: B. Fernando, A. Habrard, M. Sebban, and T. Tuytelaars, "Unsupervised
% visual domain adaptation using subspace alignment," in ICCV, 2013
% 
% Copyright 2016 Ke YAN, Tsinghua Univ. http://yanke23.com , xjed09@gmail.com

%% default parameters
pcaCoef = 0; % see ftProc_pca_tr

defParam

%% sort samples
ftSrc = ft(maSrc,:);
ftTar = ft(~maSrc,:);

%% compute and project the samples
[ftSrcNew,pcaModelS] = ftProc_pca_tr(ftSrc,[],struct('pcaCoef',pcaCoef));
[ftTarNew,pcaModelT] = ftProc_pca_tr(ftTar,[],struct('pcaCoef',pcaCoef));
d = min(size(pcaModelS.W_prj,2),size(pcaModelT.W_prj,2));
W_prjS = pcaModelS.W_prj(:,1:d);
W_prjT = pcaModelT.W_prj(:,1:d);

ftAllNew = zeros(size(ft,1),d);
ftAllNew(maSrc,:) = ftSrcNew*W_prjS'*W_prjT;
ftAllNew(~maSrc,:) = ftTarNew;

transMdl.WS = W_prjS*W_prjS'*W_prjT;
transMdl.WT = W_prjT;
transMdl.muS = pcaModelS.mu;
transMdl.muT = pcaModelT.mu;

end
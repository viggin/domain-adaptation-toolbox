function [ftAllNew,transMdl] = ftTrans_itl(ft,maSrc,target,maLabeled,param)
%Wrapper of the Information-Theoretical Learning (ITL) algorithm
% 
%	A feaure-level transfer learning (domain adaptation) algorithm which
% learns a domain-invariant subspace. Application scope:
%	+ two discrete domains (source and target)
%	+ labeled source domain
%	+ unlabeled target domain
%	+ label type: classification
%
% ftAll:	All samples in all domains. n-by-m matrix, n is the number of 
%	samples, m is the dimension of features.
% maSrc:	n-by-1 logical vector, maSrc(i)=true if i is from the source
%	domain, 0 if target domain.
% target:	The class labels of the source domain, nsrc-by-1 matrix, nsrc 
%	is the number of source	samples. target(i)=j if the i'th sample is from
%	the j'th class.
% maLabeled:	Mask for the labeled samples. n-by-1 matrix, must be the
% same with maSrc.
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
%	The ftProc_pca_tr function is a part of the PRTools toolbox.
%	Notice from the author of the paper: for practice use, you need to 
% properly preprocess the data and tune lambda and d.
% 
% ref: Y. Shi and F. Sha, "Information-theoretical learning of discriminative
% clusters for unsupervised domain adaptation," in ICML, 2012.
% 
% Copyright 2016 Ke YAN, Tsinghua Univ. http://yanke23.com , xjed09@gmail.com

addpath infometric_0.1

%% default parameters
pcaCoef = 0; % see ftProc_pca_tr
lambda = 1; % regularization parameter

defParam

%% sort samples
if any(maSrc~=maLabeled)
	error('maLabeled must be the same with maSrc')
end

ftSrc = ft(maSrc,:);
ftTar = ft(~maSrc,:);

%% compute
[~,pcaModelS] = ftProc_pca_tr(ftSrc,[],struct('pcaCoef',pcaCoef));
[~,pcaModelT] = ftProc_pca_tr(ftTar,[],struct('pcaCoef',pcaCoef));
d = min(size(pcaModelS.W_prj,2),size(pcaModelT.W_prj,2));
W_prjT = pcaModelT.W_prj(:,1:d);

L = infometric(W_prjT, ftSrc, target, ftTar, lambda);

%% project the samples
ftAllNew = zeros(size(ft,1),d);
ftAllNew(maSrc,:) = ftSrc*L;
ftAllNew(~maSrc,:) = ftTar*L;
transMdl.L = L;

end
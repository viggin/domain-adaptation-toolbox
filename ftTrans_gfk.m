function [ftAllNew,transMdl] = ftTrans_gfk(ftAll,maSrc,target,maLabeled,param)
%Wrapper of the Geodesic Flow Kernel (GFK) algorithm
% 
%	A feaure-level transfer learning (domain adaptation) algorithm which
% learns a domain-invariant subspace. Application scope:
%	+ two discrete domains (source and target)
%	+ labeled or unlabeled source domain
%	+ unlabeled target domain
%	+ label type: classification or regression
%
% ftAll:	All samples in all domains. n-by-m matrix, n is the number of 
%	samples, m is the dimension of features.
% maSrc:	n-by-1 logical vector, maSrc(i)=true if i is from the source
%	domain, 0 if target domain.
% target:	empty or nsrc-by-1 matrix, nsrc is the number of source 
%	samples. If not empty, the PLS algorithm will be used to learn
%	discriminative subspace for the source domain.
% maLabeled:	Mask for the labeled samples. n-by-1 matrix, must be the
%	same with maSrc or all 0.
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
%	The original method in ref 1 computes a kernel, but it is not very
% convenient to use, so we instead decompose the kernel to a projection 
% matrix L which is similar to ref 2.
% 
%	The GFK code can be obtained from http://www-scf.usc.edu/~boqinggo/
% ref 1: B. Gong, Y. Shi, F. Sha, and K. Grauman, "Geodesic flow kernel
% for unsupervised domain adaptation," in CVPR, 2012.
% 2: B. Gong, K. Grauman, and F. Sha, "Learning kernels for unsupervised
% domain adaptation with applications to visual object recognition", in
% IJCV, 2014
%	The plslda toolbox is needed if using PLS subspace, 
% http://www.mathworks.com/matlabcentral/fileexchange/47767-libpls-1-95-zip/
% 
%	Copyright 2016 Lu KOU, PolyU HK; Ke YAN, Tsinghua Univ. 
% http://yanke23.com , xjed09@gmail.com

addpath ToRelease_GFK

%% default parameters
dr = 0; % ratio that controls the dimension of the subspace. If 0, will be 
		% automatically computed according to ref 1. See the code
		% below.
bSup = 0; % if use the label info in source domain

defParam

%% sort samples
if any(maSrc~=maLabeled)
	error('maLabeled must be the same with maSrc')
end
ftAll(maSrc,:) = zscore(ftAll(maSrc,:)); % according to gfk's sample code
ftAll(~maSrc,:) = zscore(ftAll(~maSrc,:));

nFt = size(ftAll,2);
ftSrc = ftAll(maSrc,:);
ftTar = ftAll(~maSrc,:);

%% compute
if ~bSup
	Ps = pca(ftSrc);  % source subspace
else
	addpath plslda
	model = pls_basis(ftSrc,target,min(size(ftSrc)),'none');
	Ps = model.weight;
end
Pt = pca(ftTar);  % target subspace
Pst = pca(ftAll);

%% select subspace dimension according to ref 1
maxd = min(cellfun(@(x)size(x,2),{Ps,Pt,Pst}));
Ps = Ps(:,1:maxd); 
Pt = Pt(:,1:maxd);
Pst = Pst(:,1:maxd);
PsPst = abs(diag(Ps'*Pst)); % abs?
PtPst = abs(diag(Pt'*Pst));
PsLen = sqrt(sum(Ps.^2,1))'; % should be all 1
PtLen = sqrt(sum(Pt.^2,1))';
PstLen = sqrt(sum(Pst.^2,1))';
alphas = acos(PsPst./PsLen./PstLen);
betas = acos(PtPst./PtLen./PstLen);
Dd = (sin(alphas)+sin(betas))/2;
d = find(Dd>1-1e-3,1,'first');
if isempty(d), d = maxd; end
d = min(d,floor(nFt/2));

if dr ~= 0 % manually set the subspace dimension
	d = floor(min(maxd,floor(nFt/2))*dr);
end

%% compute
Ps = Ps(:,1:d);
Pt = Pt(:,1:d);
G = GFK([Ps,null(Ps')], Pt);
[TL,TD]=ldl(G);
L=TL*(TD.^0.5);
% 	A = chol(G+eps*20*eye(size(G,1)));
% 	L = A'; % similar to ldl
transMdl.W = real(L(:,1:2*d)); % imaginary after d*2 because rank deficient

%% project the samples
ftAllNew = ftAll*transMdl.W;

end

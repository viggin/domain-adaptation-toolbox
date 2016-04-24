function [ftAllNew,transMdl] = ftTrans_tca(ftAll,maSrc,target,maLabeled,param)
%Transfer Component Analysis (TCA)
% 
%	A feaure-level transfer learning (domain adaptation) algorithm which
% learns a domain-invariant subspace. Application scope:
%	+ two discrete domains (source and target)
%	+ labeled or unlabeled or partially labeled source domain
%	+ labeled or unlabeled or partially labeled target domain
%	+ label type: classification or regression
%
% ftAll:	All samples in all domains. n-by-m matrix, n is the number of 
%	samples, m is the dimension of features.
% maSrc:	n-by-1 logical vector, maSrc(i)=true if i is from the source
%	domain, 0 if target domain.
% target:	When some samples in any domain are labeled, their labels can
%	be provided in this variable to enhance the discriminative power of the
%	the learned features. ntr-by-1 matrix, ntr is the number of labeled 
%	samples. Classification problems should use class indices as labels, 
%	i.e. target(i)=j if the i'th labeled sample is from the j'th class.
% maLabeled:	Mask for the labeled samples. n-by-1 matrix,
%	maLabeled(i)=true if sample i is labeled, false else.
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
% ref: S. J. Pan, I. W. Tsang, J. T. Kwok, and Q. Yang, "Domain adaptation
%	via transfer component analysis," Neural Networks, IEEE Trans, 2011
% 
% Copyright 2016 Ke YAN, Tsinghua Univ. http://yanke23.com , xjed09@gmail.com

%% default parameters
isRegress = 0; % 0 for classification problem, 1 for regression
kerName = 'lin'; % kernel name, see the next cell ("kernels")
kerSigma = 10; % kernel parameter, see the next cell ("kernels")
bSstca = false; % 0 for TCA if no label information is considered for all 
	% samples, 1 for semisupervised TCA (SSTCA) if some labels are
	% considered.
doSample = false; % when there are too many unlabeled data, eigenvalue 
	% decomposition can be slow. Setting this variable to true makes the
	% code to sample some unlabeled data.
nSmpRatio = 1; % if doSample=true, the number of unlabeled data to sample 
	% will be ceil(ntr*nSmpRatio)
mu = 1; % the weight of the regularization term, see the ref
m = 30; % the dimension of the subspace

% SSTCA params
lambda = 1; % weight of the geometry term
knn = 5; % #neighbor when computing the Laplacian matrix
geoSigma = .01;
gamma = .1; % the weight of the supervised term, only useful in SMIDA, see 
	% the ref

defParam % set user-defined hyper-parameters

%% kernels
nm = @(X,p)repmat(sum(X.^2,2),1,p);
linKer = @(X1,X2)X1*X2';
rbfKer = @(X1,X2)exp(-(nm(X1,size(X2,1))+nm(X2,size(X1,1))'-2*X1*X2')/2/kerSigma^2);
polyKer = @(X1,X2)(1+kerSigma*X1*X2').^2;
lapKer = @(X1,X2)exp(-pdist2(X1,X2)/kerSigma);
if strcmpi(kerName,'lin'), kerFun = linKer;
elseif strcmpi(kerName,'poly'), kerFun = polyKer;
elseif strcmpi(kerName,'rbf'), kerFun = rbfKer;
elseif strcmpi(kerName,'lap'), kerFun = lapKer;
else error('unknown ker'); end

%% sort samples
ftLabeled = ftAll(maLabeled,:);
ftUnlabeled = ftAll(~maLabeled,:);
ftSrc = ftAll(maSrc,:);
ftTar = ftAll(~maSrc,:);
nl = size(ftLabeled,1);
nul = size(ftUnlabeled,1);
nsrc = size(ftSrc,1);
ntar = size(ftTar,1);
if doSample % sample some target domain data
	rng(0)
    ntarNew = min(ntar,floor(nsrc*nSmpRatio));
    ftAll1 = [ftSrc;ftTar(randperm(ntar,ntarNew),:)];
	ntar = ntarNew;
else
    ftAll1 = [ftSrc;ftTar];
end

%% compute
K = kerFun(ftAll1,ftAll1);
L = [ones(nsrc)/nsrc^2, -ones(nsrc,ntar)/nsrc/ntar;
    -ones(ntar,nsrc)/nsrc/ntar,ones(ntar)/ntar^2];
H = eye(nsrc+ntar)-ones(nsrc+ntar)/(nsrc+ntar);

if ~bSstca % TCA
    A = (K*L*K+mu*eye(nsrc+ntar))\K*H*K;
else % SSTCA
    % laplacian matrix
    M = squareform(pdist(ftAll1));
	geoSigma = mean(pdist(ftAll1));
	M = exp(-M.^2/2/geoSigma^2);
	M = M-eye(nsrc+ntar);
	Msort1 = sort(M,'descend');
    for p = 1:nsrc+ntar
        M(M(:,p)<Msort1(knn,p),p) = 0; % k near neighbors
    end
    M = max(M,M');
	D = diag(sum(M,2));
	Lgeo = D-M;

    if isRegress==0 && all(target==floor(target)) % clsf
        Kyy = repmat(target,1,nl)==repmat(target,1,nl)';
    else % regress
        Kyy = target*target'; % lin ker
    end
    Kyy_tilde = (1-gamma)*eye(nsrc+ntar);
	
	maLabeledNew = [maLabeled(maSrc);maLabeled(~maSrc)];
	Kyy_tilde(maLabeledNew,maLabeledNew) = ...
		Kyy_tilde(maLabeledNew,maLabeledNew)+gamma*Kyy;
	
    A = (K*(L+lambda*Lgeo)*K+mu*eye(nsrc+ntar))\...
        K*H*Kyy_tilde*H*K;
end

maBad = isnan(A)|isinf(A);
if any(any(maBad)), A(maBad)=0; end
[V,D] = eig(A);
[D,I] = sort(diag(D),'descend');
transMdl.W = real(V(:,I(1:m)));

%% project the samples
KNew = kerFun(ftAll,ftAll1);
ftAllNew = KNew*transMdl.W;

end

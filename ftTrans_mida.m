function [ftAllNew,transMdl] = ftTrans_mida(ftAll,domainFtAll,target,...
	maLabeled,param)
%Maximum Independence Domain Adaptation (MIDA)
% 
%	A feaure-level transfer learning (domain adaptation) algorithm which
% augments the features then learns a domain-invariant subspace.
% Application scope:
%	+ two or multiple discrete domains
%	+ continuous distributional change
%	+ labeled or unlabeled or partially labeled source domain
%	+ labeled or unlabeled or partially labeled target domain
%	+ label type: classification or regression
%	There is no need to distinguish which domain a sample is from. The domain
% information is contained in the variable domainFtAll, which contains the
% domain features of all samples.
%	Domain features indicate the background of samples. For example, if
% the samples are from md discrete domains, then domainFtAll can be a 
% n-by-md matrix, domainFtAll(i,j)=1 if sample i is from the j'th domain, 0
% else. If the samples have a time order and their distribution changes
% continuously along with time, then domainFtAll can be n-by-1 and
% domainFtAll(i) is the "time" of sample i. Multiple background info can be
% integrated into domainFtAll in the similar way. See the ref for details.
%
% ftAll:	All samples in all domains. n-by-m matrix, n is the number of 
%	samples, m is the dimension of features.
% domainFtAll:	Domain features of all samples. n-by-md matrix, md is the 
%	dimension of domain features.
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

% ftAllNew:	All samples in the learned subspace.
% transMdl:	A struct containing the model, transMdl.W is the projection
%	matrix.

% ref: Ke Yan, Lu Kou, and David Zhang, "Domain Adaptation via Maximum 
%	Independence of Domain Features," http://arxiv.org/abs/1603.04535
% Copyright 2016 Ke YAN, Tsinghua Univ. http://yanke23.com , xjed09@gmail.com

%% default parameters
isRegress = 0; % 0 for classification problem, 1 for regression
kerName = 'lin'; % kernel name, see the next cell ("kernels")
kerSigma = 10; % kernel parameter, see the next cell ("kernels")
bSmida = true; % 0 for MIDA if no label information is considered for all 
	% samples, 1 for semisupervised MIDA (SMIDA) if some labels are
	% considered.
doSample = false; % when there are too many unlabeled data, eigenvalue 
	% decomposition can be slow. Setting this variable to true makes the
	% code to sample some unlabeled data.
nSmpRatio = 1; % if doSample=true, the number of unlabeled data to sample 
	% will be ceil(ntr*nSmpRatio)
mu = 1; % the weight of the variance term, see the ref
m = 30; % the dimension of the subspace
gamma = .1; % the weight of the supervised term, only useful in SMIDA, see 
	% the ref
ftAugType = 1; % feature augmentation, 0: no aug; 1: aug with domainFt; 
	% 2: frustratingly easy aug, only for discrete domains, see the ref

defParam % set user-defined hyper-parameters

%% kernels
nm = @(X,p)repmat(sum(X.^2,2),1,p);
linKer = @(X1,X2)X1*X2';
rbfKer = @(X1,X2)exp(-(nm(X1,size(X2,1))+nm(X2,size(X1,1))'-2*X1*X2')/2/kerSigma^2);
lapKer = @(X1,X2)exp(-pdist2(X1,X2)/kerSigma);
polyKer = @(X1,X2)(1+kerSigma*X1*X2').^2;
if strcmpi(kerName,'lin'), kerFun = linKer; % linear kernel
elseif strcmpi(kerName,'poly'), kerFun = polyKer; % polynomial kernel
elseif strcmpi(kerName,'rbf'), kerFun = rbfKer;
elseif strcmpi(kerName,'lap'), kerFun = lapKer; % Laplacian kernel
else error('unknown kernel'); end

%% sort samples
ntr = nnz(maLabeled);
nAll = size(ftAll,1);
Df = double(domainFtAll);
% Df = zscore(domainFtAll); % bad
if doSample && nSmpRatio<inf % sample some unlabeled data
	rng(0)
	nUnlabeledUsed = min(nAll-ntr,ceil(ntr*nSmpRatio));
	id = randperm(nAll-ntr,nUnlabeledUsed); % rand sel
	idTest = find(~maLabeled);
	ftUsed = [ftAll(maLabeled,:);ftAll(idTest(id),:)];
	dfUsed = [Df(maLabeled,:);Df(idTest(id),:)];
	nUsed = ntr+nUnlabeledUsed;
else
	ftUsed = [ftAll(maLabeled,:);ftAll(~maLabeled,:)];
	dfUsed = [Df(maLabeled,:);Df(~maLabeled,:)];
	nUsed = nAll;
end
Kd = linKer(dfUsed,dfUsed);

%% feature augmentation
if ftAugType==0
	ftUsedAug = ftUsed;
elseif ftAugType==1
	ftUsedAug = [ftUsed,dfUsed];
elseif ftAugType==2
	nFt = size(ftUsed,2);
	nDomain = size(dfUsed,2);
	ftUsedAug = zeros(nUsed,nFt*(nDomain+1));
	ftUsedAug(:,1:nFt) = ftUsed;
	for p = 1:nDomain
		ftUsedAug(dfUsed(:,p)==1,nFt*p+1:nFt*(p+1)) = ftUsed(dfUsed(:,p)==1,:);
	end
end

%% compute
Kx = kerFun(ftUsedAug,ftUsedAug);
H = eye(nUsed)-ones(nUsed)/(nUsed);

if ~bSmida % MIDA
	
	A = Kx*(mu*H-H*Kd*H)*Kx;
	
else % SMIDA
	
	target = double(target);
	if isRegress==0 && all(target==floor(target)) % classification
		Ky = repmat(target,1,ntr)==repmat(target,1,ntr)';
% 		Ky = (Kyy-.5)*2; % no need
	else % regression
		target = zscore(target);
		Ky = target*target'; % linear ker
	end
	
	Kyy_tilde = zeros(nUsed);
	Kyy_tilde(1:ntr,1:ntr) = Ky;
	
	A = Kx*(H*(-Kd+gamma*Kyy_tilde)*H+mu*H)*Kx;
end

A = (A+A')/2; % to compensate float num error
[V,D] = eig(A);
[D,I] = sort(diag(real(D)),'descend');
transMdl.W = real(V(:,I(1:m)));

%% project the samples
if ftAugType==0
	ftAllAug = ftAll;
elseif ftAugType==1
	ftAllAug = [ftAll,Df];
elseif ftAugType==2
	ftAllAug = zeros(nUsed,nFt*(nDomain+1));
	ftAllAug(:,1:nFt) = ftAll;
	for p = 1:nDomain
		ftAllAug(Df(:,p)==1,nFt*p+1:nFt*(p+1)) = ftAll(Df(:,p)==1,:);
	end
end

KxAll = kerFun(ftAllAug,ftUsedAug);
ftAllNew = KxAll*transMdl.W;

end

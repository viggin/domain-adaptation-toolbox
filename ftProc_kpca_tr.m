function [Xnew, model] = ftProc_kpca_tr(X,Y,param)
%ftProc_kpca_tr Training Kernal Principal Component Analysis (KPCA).
%	[Xnew, model] = ftProc_kpca_tr(X,Y,param)
%	X: a matrix, each sample is a row. Y is useless.
% 	XNEW is a matrix, each row is the PCA feature of a sample in X.
%	MODEL is the struct containing coefficients of PCA.
%	PARAM is a struct of the parameters for PCA, it can contain:
%		pcaCoef: if 0, extract all principal components;
%			elseif integer, extract pcaCoef PCs;
%			else >0 and <1, extract PC with energy ratio of pcaCoef. Default 0.
%		kernel:	'lin','poly','rbf','lap'. Default 'poly'.
%		gamma and degree: kernel parameters. See the code.

%	Copyright 2015 Ke YAN, Tsinghua Univ. http://yanke23.tk, xjed09@gmail.com

pcaCoef = 0;
kernel = 'poly';
gamma = 1;
degree = 2;

defParam

nm = @(X,p)repmat(sum(X.^2,2),1,p);
linKer = @(X1,X2)X1*X2';
rbfKer = @(X1,X2)exp(-(nm(X1,size(X2,1))+nm(X2,size(X1,1))'-2*X1*X2')/2/gamma^2);
lapKer = @(X1,X2)exp(-pdist2(X1,X2)/gamma);
polyKer = @(X1,X2)(1+gamma*X1*X2').^degree;
if strcmpi(kernel,'lin'), kerFun = linKer;
elseif strcmpi(kernel,'poly'), kerFun = polyKer;
elseif strcmpi(kernel,'rbf'), kerFun = rbfKer;
elseif strcmpi(kernel,'lap'), kerFun = lapKer; % Laplacian
else error('unknown ker'); end

K = kerFun(X,X);
[nSmp,nFt] = size(X);
H = ones(nSmp)/nSmp;
Kcentered = K-K*H-H*K+H*K*H;
[V,D] = eig(Kcentered);

D = real(diag(D));
[evs,I] = sort(D,'descend');
V = real(V(:,I));

D1 = cumsum(evs);
D1 = D1/D1(end);
% calculate the dimension of the data after PCA projection
if pcaCoef == 0
	postPcaDim = min(nSmp,nFt);
elseif pcaCoef > 0 && pcaCoef < 1
	pcaRate = pcaCoef;
	postPcaDim = nnz(D1 < pcaRate)+1;
elseif pcaCoef >= 1 && floor(pcaCoef) == pcaCoef
	postPcaDim = pcaCoef;
end

% normalize the projection vectors with nonzero eigenvalues
W_prj = V(:,1:postPcaDim);
nzev = evs(1:postPcaDim)>1e-6;
W_prj(:,nzev) = W_prj(:,nzev)./repmat(sqrt(evs(nzev))',nSmp,1);

% project X to KPCA subspace
Xnew = K*W_prj;

model.Xtrain = X;
model.W_prj = W_prj;
model.eigVals = evs;
model.postPcaDim = postPcaDim;
model.kerFun = kerFun;

end

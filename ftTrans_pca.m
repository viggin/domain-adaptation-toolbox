function [ftAllNew,transMdl] = ftTrans_pca(ftAll,maSrc,target,maLabeled,param)
%A baseline to other domain adaptation methods, use PCA or KPCA to extract
% features.
% 
% Application scope:
%	+ Any scenario since no domain adaptation is done.
% 
% ftAll:	All samples in all domains. n-by-m matrix, n is the number of 
%	samples, m is the dimension of features.
% maSrc,target,maLabeled:	useless.
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
% The ftProc_(k)pca_tr function is a part of the PRTools toolbox.
% if kerName is set in param, use KPCA; else, use original linear PCA.
% 
% Copyright 2016 Ke YAN, Tsinghua Univ. http://yanke23.com , xjed09@gmail.com

%% default parameters
pcaCoef = 0; % see ftProc_pca_tr
kerName = 'poly'; % see ftProc_kpca_tr
kerSigma = 1;
degree = 2;

defParam

%% compute and project the samples
if ~isfield(param,'kerName')
	[ftAllNew,pcaModel] = ftProc_pca_tr(ftAll,[],struct('pcaCoef',pcaCoef));
else
	[ftAllNew,pcaModel] = ftProc_kpca_tr(ftAll,[],struct('pcaCoef',pcaCoef,...
		'kernel',kerName,'kerSigma',kerSigma,'degree',degree));
end

transMdl = pcaModel;

end
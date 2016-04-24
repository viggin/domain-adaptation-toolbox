function [ftAllNew,transMdl] = ftTrans_ssa(ftAll,maSrc,target,maLabeled,param)
%Wrapper of the Stationary Subspace Analysis (SSA) algorithm
% 
%	A feaure-level transfer learning (domain adaptation) algorithm which
% learns a stationary subspace. Application scope:
%	+ continuous distributional change
%	+ unlabeled samples
%	Samples in ftAll are assumed to be a multivariate time series.
%
% ftAll:	All samples. n-by-m matrix, n is the number of 
%	samples, m is the dimension of features.
% maSrc,target,maLabeled: useless
% 
% param:	Struct of hyper-parameters, please see the first cell of this
%	program ("default parameters") for details. You can set parameter p to 
%	x by setting param.p = x. For parameters that are not set, default 
%	values will be used.

% ftAllNew:	All samples in the learned subspace.
% transMdl:	A struct containing the model, transMdl.W is the projection
%	matrix.
% 
% Please download the ssa toolbox at: http://mloss.org/revision/download/851/
% ref: P. Von Bunau, et al, "Finding stationary subspaces in multivariate 
% time series," Physical review letters, 2009.
% Copyright 2016 Ke YAN, Tsinghua Univ. http://yanke23.com , xjed09@gmail.com

addpath ssa_toolbox-1.3

%% default parameters
d = 0; % num of features to return, if 0, will be set to m
m = 1; % num of stationary components
% other params see the doc of the function ssa

defParam

%% compute
if d==0, d = m; end
[est_Ps, est_Pn, est_As, est_An, ssa_results] = ssa(ftAll', m);

%% project the samples
ftAllNew = [ssa_results.s_src;ssa_results.n_src]';
ftAllNew = ftAllNew(:,1:d);
transMdl = ssa_results;

end

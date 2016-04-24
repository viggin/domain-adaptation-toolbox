function [pred,model,prob] = mdlTrans_lapsvm(Xl,Yl,Xu,param)
%Wrapper of Laplacian SVM
% 
%	This is a model-level semi-supervised (transductive) learning
% algorithm. It may be used as a baseline to other domain adaptation 
% methods.
% Application scope:
%	+ partially labeled data
%	+ label type: classification

% Xl:	labeled X
% Yl:	labels of Xl, should be a vector containing class indices starting from 1
% Xu:	unlabeled X
% 
% param: Struct of hyper-parameters, please see the first cell of this
%	program ("default parameters") for details. You can set parameter p to 
%	x by setting param.p = x. For parameters that are not set, default 
%	values will be used.

% pred : predicted labels for Xu
% prob : the confidence of pred
% model: model. For multi-class problems, one-vs-all strategy will be used
%	and model{i} is the model for class i
%
%	The lapsvm toolbox can be obtained from 
% https://github.com/tknandu/LapTwinSVM/tree/master/Primal_LapSVM/lapsvmp_v02
% ref: M. Belkin, P. Niyogi, and V. Sindhwani, "Manifold regularization: A
% geometric framework for learning from labeled and unlabeled examples,"
% J. Mach. Learn. Res., 2006.
% 
% Copyright 2016 Ke YAN, Tsinghua Univ. http://yanke23.com , xjed09@gmail.com

addpath(genpath('lapsvmp_v02'))

%% default parameters
% please find the details and other parameters in functions lapsvmp, 
% calckernel, and laplacian
t = 2;
% c = 1;
g = .001;
gamma_I = 1;
gamma_A = 1e-5;
knn = 5;

defParam

%% kernels
if t == 0, Kernel = 'linear'; KernelParam = 0;
elseif t == 1, Kernel = 'poly'; KernelParam = g;
elseif t == 2, Kernel = 'rbf'; KernelParam = sqrt(1/2/g);
else error('wrong ker'); end

%% make options
options = make_options('gamma_I',gamma_I,'gamma_A',gamma_A,'NN',knn,...
	'Kernel',Kernel,'KernelParam',KernelParam);
options.Verbose = 0;
options.UseBias = 1;
options.UseHinge = 1;
options.LaplacianNormalize = 0;
options.NewtonLineSearch = 0;
if t == 1,options.polyDegree = 2;end

%% create the 'data' structure
X = [Xl;Xu];
data.X = X;
data.K = calckernel(options,X,X);
data.L = laplacian(options,X);

%% train the classifier and predict
options.Cg = 1; % PCG
options.MaxIter = 1000; % upper bound
options.CgStopType = 1; % 'stability' early stop
options.CgStopParam = 0.015; % tolerance: 1.5%
options.CgStopIter = 3; % check stability every 3 iterations

nCls = max(Yl);
nTest = size(Xu,1);
if nCls > 2 % multi-class
	warning('currently using one-vs-all strategy on lapSVM, one-vs-one may be better.');
	prob = nan(nTest,nCls);
	for iCls = 1:nCls
		Y1 = Yl == iCls;
		data.Y = [(Y1-.5)*2;zeros(nTest,1)];% {0,1} to {-1,0,+1}
		classifier = lapsvmp(options,data);
		% fprintf('It took %f seconds.\n',classifier.traintime);
		
		out = data.K(:,classifier.svs)*classifier.alpha+classifier.b;
		prob(:,iCls) = out(end-nTest+1:end);
		model{iCls} = classifier;
	end
	[~,pred] = max(prob,[],2);
	
else
	data.Y = [(Yl-1.5)*2;zeros(nTest,1)];
	classifier = lapsvmp(options,data);
	
	out = data.K(:,classifier.svs)*classifier.alpha+classifier.b;
	prob = out(end-nTest+1:end);
	pred = sign(prob)/2+1.5;
	model = classifier;
end

end
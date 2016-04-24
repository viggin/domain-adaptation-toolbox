function [pred,model] = mdlTrans_lapridge(Xl,Yl,Xu,param)
%Laplacian ridge regression
% 
%	This is a model-level semi-supervised (transductive) learning
% algorithm. It may be used as a baseline to other domain adaptation 
% methods.
% Application scope:
%	+ partially labeled data
%	+ label type: regression

% Xl:	labeled X
% Yl:	labels of X
% Xu:	unlabeled X
% 
% param: Struct of hyper-parameters, please see the first cell of this
%	program ("default parameters") for details. You can set parameter p to 
%	x by setting param.p = x. For parameters that are not set, default 
%	values will be used.
% 
% pred : predicted labels for Xu
% 
% ref: M. Belkin, P. Niyogi, and V. Sindhwani, "Manifold regularization: A
% geometric framework for learning from labeled and unlabeled examples,"
% J. Mach. Learn. Res., 2006.
% Copyright 2016 Ke YAN, Tsinghua Univ. http://yanke23.com , xjed09@gmail.com

addpath(genpath('lapsvmp_v02'))

%% default parameters
% please find the details and other parameters in functions calckernel and laplacian
t = 2;
g = .001;
gamma_I = 1;
gamma_A = 1e-5;
knn = 5;

defParam

%% kernels
if t==0, Kernel = 'linear'; KernelParam = 0;
elseif t==1, Kernel = 'poly'; KernelParam = g;
elseif t==2, Kernel = 'rbf'; KernelParam = sqrt(1/2/g);
else error('wrong ker'); end

%% make options
options = make_options('gamma_I',gamma_I,'gamma_A',gamma_A,'NN',knn,...
	'Kernel',Kernel,'KernelParam',KernelParam);
options.LaplacianNormalize = 0;
if t==1,options.polyDegree = 2;end

%% sort data
nTrain = size(Xl,1);
nTest = size(Xu,1);
Xl = [ones(nTrain,1),Xl];
Xu = [ones(nTest,1),Xu];

%% compute
X = [Xl;Xu];
K = calckernel(options,X,X);
L = laplacian(options,X);

J = diag([ones(1,nTrain),zeros(1,nTest)]);
I = eye(nTrain+nTest);
Yl = [Yl;zeros(nTest,1)];
alpha = (J*K + gamma_A*nTrain*I + gamma_I*nTrain/(nTrain+nTest)^2*L*K) \ Yl;

%% predict
out = K*alpha;
pred = out(nTrain+1:nTrain+nTest);
model.K = K;
model.alpha = alpha;

end

function [pred,model,prob] = mdlTrans_tsvm(Xl,Yl,Xu,param)
%Wrapper of Transducive SVM
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
%	The svm-light toolbox and matlab interface can be obtained from
% http://svmlight.joachims.org/ and https://github.com/sods/svml
% ref: T. Joachims, "Transductive inference for text classification using 
% support vector machines," 1999
% 
% Copyright 2016 Ke YAN, Tsinghua Univ. http://yanke23.com , xjed09@gmail.com

svmlPath = 'svml-master';
addpath(svmlPath)

%% default parameters
% please find the details and other parameters in function svml
t = 2;
c = 1;
g = .001;
TransPosFrac = nan; % Fraction of unlabeled examples to be classified into 
					% the positive class. If nan, we will estimate it using
					% Xu
defParam

%% train the classifier and predict
X = [Xl;Xu];
nTest = size(Xu,1);

nCls = max(Yl);
if nCls > 2
	warning('currently using one-vs-all strategy on TSVM, one-vs-one may be better.');
	prob = nan(nTest,nCls);

	for iCls = 1:nCls
		Y1 = Yl==iCls;
		Y1 = (Y1-.5)*2;% {0,1} to {-1,+1}
		if isnan(TransPosFrac)
			TransPosFrac = nnz(Yl==iCls)/length(Yl);
		end

		net = svml('', 'Kernel', t, 'C', c, 'KernelParam',g,'TransPosFrac',TransPosFrac,...
			'ExecPath',svmlPath,'Verbosity',0);
		Y1 = [Y1;zeros(nTest,1)]; % label 0 for the unlabeled samples
		net = svmltrain(net, X, Y1);
		prob(:,iCls) = svmlfwd(net, Xu);
		model{iCls} = net;
	end
	[~,pred] = max(prob,[],2);
	
else
	Y1 = (Yl-1.5)*2;
	if isnan(TransPosFrac)
		TransPosFrac = nnz(Y1==1)/length(Y1);
	end

	net = svml('', 'Kernel', t, 'C', c, 'KernelParam',g,'TransPosFrac',TransPosFrac,...
		'ExecPath',svmlPath,'Verbosity',0);
	Y1 = [Y1;zeros(nTest,1)];
	net = svmltrain(net, X, Y1);
	prob = svmlfwd(net, Xu);
	pred = sign(prob)/2+1.5;
	model = net;
end

end

%   Accepted options are:
%   Field      SVM light option  Range, description
%   'Verbosity'      -v       {0 .. 3}, default value 1
%                             Verbosity level
%   'Regression'     -z       {0, 1}, default value 0
%                             Switch between regression [1] and
%                             classification [0]
%   'C'              -c       (0, Inf), default value (avg. x*x)^-1
%                             Trade-off between error and margin
%   'TubeWidth'      -w       (0, Inf), default value 0.1
%                             Epsilon width of tube for regression
%   'CostFactor'     -j       (0, Inf), default value 1
%                             Cost-Factor by which training errors on
%                             positive examples outweight errors on
%                             negative examples
%   'Biased'         -b       {0, 1}, default value 1
%                             Use biased hyperplane x*w+b0 [1] instead of
%                             unbiased x*w0 [0]
%   'RemoveIncons'   -i       {0, 1}, default value 0
%                             Remove inconsistent training examples and
%                             retrain
%   'ComputeLOO'     -x       {0, 1}, default value 0
%                             Compute leave-one-out estimates [1]
%   'XialphaRho'     -o       )0, 2), default value 1.0
%                             Value of rho for XiAlpha-estimator and for
%                             pruning leave-one-out computation
%   'XialphaDepth'   -k       {0..100}, default value 0
%                             Search depth for extended XiAlpha-estimator 
%   'TransPosFrac'   -p       (0..1), default value ratio of
%                             positive and negative examples in the
%                             training data. Fraction of unlabeled
%                             examples to be classified into the positive
%                             class
%   'Kernel'         -t       {0..4}, default value 1
%                             Type of kernel function:
%                             0: linear
%                             1: polynomial (s a*b+c)^d
%                             2: radial basis function exp(-gamma ||a-b||^2)
%                             3: sigmoid tanh(s a*b + c)
%                             4: user defined kernel from kernel.h
%   'KernelParam'    -d, -g, -s, -r, -u
%                             Depending on the kernel, this vector
%                             contains [d] for polynomial kernel, [gamma]
%                             for RBF, [s, c] for tanh kernel, string for
%                             user-defined kernel
%   'MaximumQP'      -q       {2..}, default value 10
%                             Maximum size of QP-subproblems
%   'NewVariables'   -n       {2..}, default value is the value chosen
%                             for 'MaximumQP'. Number of new variables
%                             entering the working set in each
%                             iteration. Use smaller values to prevent
%                             zig-zagging
%   'CacheSize'      -m       (5..Inf), default value 40.
%                             Size of cache for kernel evaluations in MB
%   'EpsTermin'      -e       (0..Inf), default value 0.001
%                             Allow that error for termination criterion
%                             [y [w*x+b] - 1] < eps
%   'ShrinkIter'     -h       {5..Inf}, default value 100.
%                             Number of iterations a variable needs to be
%                             optimal before considered for shrinking
%   'ShrinkCheck'    -f       {0, 1}, default value 1
%                             Do final optimality check for variables
%                             removed by shrinking. Although this test is
%                             usually positive, there is no guarantee
%                             that the optimum was found if the test is
%                             omitted.
%   'TransLabelFile' -l       String. File to write predicted labels of
%                             unlabeled examples into after transductive
%                             learning.
%   'AlphaFile'      -a       String. Write all alphas to this file after
%                             learning (in the same order as in the
%                             training set).

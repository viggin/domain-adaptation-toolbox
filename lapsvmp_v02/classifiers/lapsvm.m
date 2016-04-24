function classifier = lapsvm(options,data)
% {lapsvm} trains a Laplacian SVM classifier (dual solved with Libsvm).
%     
%      classifier = lapsvm(options,data)
%
%      options: a structure with the following fields
%               options.gamma_A: regularization parameter (ambient space)
%               options.gamma_I: regularization parameter (intrinsic norm)
%               
%               [optional fields]
%               options.UseBias: {0,1} i.e. use or not a bias (fx=Kalpha+b)               
%      data: a structure with the following fields
%            data.X: a M-by-D matrix of M D-dimensional training examples
%            data.K: a M-by-M kernel Gram matrix of M training examples
%            data.Y: a M-by-1 label vector in {-1,0,+1}, where 0=unlabeled
%            data.L: a M-by-M matrix of the graph Laplacian            
%
%      classifier: structure of the trained classifier (see the
%                  'saveclassfier' function)
%
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it
%         * based on the code of Vikas Sindhwani, vikas.sindhwani@gmail.com

tic
if ~isfield(options,'UseBias'),           options.UseBias=1; end

C=1;
parameters =    [4 ... % kernel_type (4=Gram Matrix)
                1 ... % deg
                1 ... % gamma
                0 ... % coef0
                C ... % C
                40.00 ... % cache
                0.001 ... % eps
                0 ... % svm_type
                0.5 ... % nu
                0.1 ... % p
                1 ... % shrinking
                ];

I=eye(size(data.K,1));
lab=(data.Y~=0);
l=nnz(lab);

if isempty(data.L) || options.gamma_I==0
    G=I/(2*options.gamma_A);
else
    G=(2*options.gamma_A*I + 2*options.gamma_I*data.L*data.K)\I;
end

Gram=data.K(lab,:)*G(:,lab);
Ylab=data.Y(lab);

[beta, svs, b, nsv, nlab] = mexGramSVMTrain(Gram', Ylab', parameters);

if nlab(1)==-1
    beta=-beta;
else
    b=-b;
end

betaz=zeros(l,1);
betaz(svs)=beta';
alpha=G(:,lab)*betaz;
sec=toc;

classifier = saveclassifier('lapsvm',1:length(data.Y),alpha, ...
                            data.X,b*options.UseBias,options,sec);
                        
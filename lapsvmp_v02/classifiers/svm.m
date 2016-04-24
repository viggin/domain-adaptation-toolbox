function classifier = svm(options,data)
% {svm} trains a SVM classifier (using Libsvm).
%     
%      classifier = svm(options,data)
%
%      options: a structure with the following fields
%               options.gamma_A: regularization parameter (ambient space)
%
%               [optional fields]
%               options.UseBias: {0,1} i.e. use or not a bias.              
%      data: a structure with the following fields
%            data.X: a N-by-D matrix of N D-dimensional training examples
%            data.K: a N-by-N kernel Gram matrix of N training examples
%            data.Y: a N-by-1 label vector in {-1,+1}
%
%      classifier: structure of the trained classifier (see the
%                  'saveclassfier' function)
%
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it
%         * based on the code of Vikas Sindhwani, vikas.sindhwani@gmail.com

tic
if ~isfield(options,'UseBias'),           options.UseBias=1; end

C=1/(2*options.gamma_A);
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

[alpha, svs, b, ~, nlab] = mexGramSVMTrain(data.K',data.Y',parameters);

alpha=alpha'; % libsvm does some weird label switching
if nlab(1)==-1 
    alpha=-alpha;
else
    b=-b;
end
sec=toc;

classifier = saveclassifier('svm',svs,alpha,data.X(svs,:),...
                            b*options.UseBias,options,sec);

% setting default paths
setpaths

% generating a random dataset
fprintf('Generating 1000 random data points...\n');
X=[randn(500,5); 2*randn(500,5)+2];
Y=[ones(500,1); -ones(500,1)];

% generating default options
options=make_options('gamma_I',1,'gamma_A',1e-5,'NN',6,'KernelParam',0.35);
options.Verbose=1;
options.UseBias=1;
options.UseHinge=1;
options.LaplacianNormalize=0;
options.NewtonLineSearch=0;

% creating the 'data' structure
data.X=X;
data.Y=zeros(size(Y));
data.Y(1:50)=1; % 50 labeled points of class +1
data.Y(501:550)=-1; % 50 labeled points of class -1

fprintf('Computing Gram matrix and Laplacian...\n\n');
data.K=calckernel(options,X,X);
data.L=laplacian(options,X);

% training the classifier
fprintf('Training LapSVM in the primal with Newton''s method...\n');
classifier=lapsvmp(options,data);

% computing error rate
fprintf('It took %f seconds.\n',classifier.traintime);
out=sign(data.K(:,classifier.svs)*classifier.alpha+classifier.b);
er=100*(length(data.Y)-nnz(out==Y))/length(data.Y);
fprintf('Error rate=%.1f\n\n',er);

% training the classifier
fprintf('Training LapSVM in the primal with early stopped PCG...\n');
options.Cg=1; % PCG
options.MaxIter=1000; % upper bound
options.CgStopType=1; % 'stability' early stop
options.CgStopParam=0.015; % tolerance: 1.5%
options.CgStopIter=3; % check stability every 3 iterations
classifier=lapsvmp(options,data);
fprintf('It took %f seconds.\n',classifier.traintime);

% computing error rate
out=sign(data.K(:,classifier.svs)*classifier.alpha+classifier.b);
er=100*(length(data.Y)-nnz(out==Y))/length(data.Y);
fprintf('Error rate=%.1f\n\n',er);


function Ex_GFK()
% This shows how to use GFK in a 1-nearest neighbor classifier.

% ref: Geodesic Flow Kernel for Unsupervised Domain Adaptation.  
% B. Gong, Y. Shi, F. Sha, and K. Grauman.  
% Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR), Providence, RI, June 2012.

% Contact: Boqing Gong (boqinggo@usc.edu)


%-------------------------I. setup source/target domains----------------------
% Four domains: { Caltech10, amazon, webcam, dslr }
src = 'webcam';
tgt = 'Caltech10';

d = 20; % subspace dimension, the following dims are used in the paper:
% webcam-dslr: 10
% dslr-amazon: 20
% webcam-amazon: 10
% caltech-webcam: 20
% caltech-dslr: 10
% caltech-amazon: 20
% Note the dim from X to Y is the same as that from Y to X.

nPerClass = 20; 
% 20 per class when Caltech/Amazon/Webcam is the source domain, and 
% 8 when DSLR is the source domain.

%--------------------II. prepare data--------------------------------------
load(['data/' src '_SURF_L10.mat']);     % source domain
fts = fts ./ repmat(sum(fts,2),1,size(fts,2)); 
Xs = zscore(fts,1);    clear fts
Ys = labels;           clear labels
Ps = princomp(Xs);  % source subspace

load(['data/' tgt '_SURF_L10.mat']);     % target domain
fts = fts ./ repmat(sum(fts,2),1,size(fts,2)); 
Xt = zscore(fts,1);     clear fts
Yt = labels;            clear labels
% Pt = princomp(Xt);  % target subspace
Pt = pca(Xt);  % target subspace
% Pt = ftProc_pca_tr(Xt);  % target subspace

fprintf('\nsource (%s) --> target (%s):\n', src, tgt);
fprintf('round     accuracy\n');
%--------------------III. run experiments----------------------------------
round = 20; % 20 random trials
tot = 0;
for iter = 1 : round 
    fprintf('%4d', iter);
    
    inds = split(Ys, nPerClass);
    Xr = Xs(inds,:);
    Yr = Ys(inds);

    %---------------III.A. PLS --------------------------------------------
    % Ps = PLS(Xr, OneOfKEncoding(Yr), 3*d);   
    % PLS generally leads to better performance.
    % A nice implementation is publicaly available at http://www.utd.edu/~herve/
    
    G = GFK([Ps,null(Ps')], Pt(:,1:d));
    [~, accy] = my_kernel_knn(G, Xr, Yr, Xt, Yt);   
    fprintf('\t\t%2.2f%%\n', accy*100);
    tot = tot + accy;
end
fprintf('mean accuracy: %2.2f%%\n\n', tot/round*100);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prediction accuracy] = my_kernel_knn(M, Xr, Yr, Xt, Yt)
dist = repmat(diag(Xr*M*Xr'),1,length(Yt)) ...
    + repmat(diag(Xt*M*Xt')',length(Yr),1)...
    - 2*Xr*M*Xt';
[~, minIDX] = min(dist);
prediction = Yr(minIDX);
accuracy = sum( prediction==Yt ) / length(Yt); 


function [idx1 idx2] = split(Y,nPerClass, ratio)
% [idx1 idx2] = split(X,Y,nPerClass)
idx1 = [];  idx2 = [];
for C = 1 : max(Y)
    idx = find(Y == C);
    rn = randperm(length(idx));
    if exist('ratio')
        nPerClass = floor(length(idx)*ratio);
    end
    idx1 = [idx1; idx( rn(1:min(nPerClass,length(idx))) ) ];
    idx2 = [idx2; idx( rn(min(nPerClass,length(idx))+1:end) ) ];
end

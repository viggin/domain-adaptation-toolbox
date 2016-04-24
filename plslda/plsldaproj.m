function [ScoresTest]=plsldaproj(model,Xtest,flag)
%+++ To project the future samples onto the training sample space.
%+++ pls_lda in this software package.
%+++ Nov.16,2007. H.D. Li.

if nargin<3;flag=0;end
[M,N]=size(Xtest);
%scale
method=model.method;
weight=model.weight;
sscale=model.scale_para;
for i=1:N; Xtest(:,i)=(Xtest(:,i)-sscale(1,i))/sscale(2,i);end
ScoresTest=Xtest*weight;
ScoresTrain=model.Xscores;
% X=[ScoresTrain;ScoresTest];
% y=[model.yreal;2*ones(M,1)];
% classplot2(X,y,flag);
% xlabel('PLS-1');
% ylabel('PLS-2');


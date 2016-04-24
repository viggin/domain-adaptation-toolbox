function [ypred,error]=plsldaval(model,Xtest,ytest)
%+++ To predict the future sample use the PLS model obtained by the function
%+++ pls_lda in this software package.
%+++ Nov.16,2007. H.D. Li.

[M,N]=size(Xtest);
%scale
method=model.method;
sscale=model.scale_para;
for i=1:N; Xtest(:,i)=(Xtest(:,i)-sscale(1,i))/sscale(2,i);end
%%%%%%%%%% predict %%%%%%%%%%%%%%%%%%%
coef=model.coef_lda_origin;
ypred=Xtest*coef(1:end-1)+coef(end);
error=sum(sign(ypred)~=ytest)/length(ytest);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





















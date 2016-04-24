function F=roccurve(yreal,ypred,flag)
%+++ yreal: with elements 1 or -1;
%+++ ypred: real values.
%+++ flag: 1: plot
%          0: no plot.
%+++ June 11,2008.

if nargin<3;flag=0;end;

yreal=sign(yreal);
thitamin=min(ypred);thitamax=max(ypred);
K=128;
thita=linspace(thitamin,thitamax,K);
Result=zeros(K,2);
for i=1:K
  r=sesp(yreal,ypred-thita(i));
  Result(i,:)=[1-r.specificity r.sensitivity];    
end
auc=abs(trapz(Result(:,1),Result(:,2)));
if flag==1
  plot(Result(:,1),Result(:,2));
  xlabel('1-specificity');ylabel('sensitivity');
end
r=sesp(yreal,ypred);
%+++ OUTPUT
F.value=Result;
F.sensitivity=r.sensitivity;
F.specificity=r.specificity;
F.accuracy=r.accuracy;
F.AUC=auc; % area under curve.


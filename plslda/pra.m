function PRA=pra(X,y,ratio,method,Nmcs,OPT)
%+++ Personalized risk assessment based on linear discriminant analysis
%+++ Input:  X: m x n  (Sample matrix)
%            y: m x 1  (measured property)
%       method: pretreatment method. Contains: autoscaling, center etc.
%         Nmcs: The number of Monte Carlo sampling.
%          OPT: =1 Print process.
%               =0 No print.
%               pareto,minmax,center or none.
%+++ Output: Structural data: CV
%+++ Hongdong Li, Jul. 4, 2010.
%+++ Revised on Dec.3, 2009.


if nargin<5;OPT=1;end;
if nargin<4;Nmcs=100;end;
if nargin<3;method='center';end;



[Mx,Nx]=size(X);
FIT=nan(Nmcs,Mx);
PRED=nan(Nmcs,Mx);
Q=ceil(Mx*ratio);

%+++ Main loop
for i=1:Nmcs
    
  perm=randperm(Mx);
  calk=perm(1:Q);
  testk=perm(Q+1:end);
  Xcal=X(calk,:);ycal=y(calk);
  Xtest=X(testk,:);ytest=y(testk); 
  
 
  C=ldapinv(Xcal,ycal,0);
  yfit=Xcal*C(1:end-1)+C(end);
  ypred=Xtest*C(1:end-1)+C(end);
  FIT(i,calk)=yfit;
  PRED(i,testk)=ypred;
  fprintf('The %dth MCS finished.\n',i);
end

%+++ output
PRA.FIT=FIT;
PRA.PRED=PRED;




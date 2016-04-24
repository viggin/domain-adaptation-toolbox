function DCV=plsldardcv(X,y,A,K,method,Nmcs,OPT,order)
%+++ Repeaded double cross validation Cross-validation for PLS-LDA
%+++ Input:  X: m x n  (Sample matrix)
%            y: m x 1  (measured property)
%            A: The max PC for cross-validation
%            K: fold. when K = m, it is leave-one-out CV
%       method: pretreatment method. Contains: autoscaling, center etc.
%          OPT: =1 Print process.
%               =0 No print.
%               pareto,minmax,center or none.
%+++ Order: =1  sorted. For CV partition.
%           =0  random, default. 
%+++ Output: Structural data: CV
%+++ Hongdong Li, Oct. 16, 2008.
%+++ Revised in Jan.12, 2009.


if nargin<8;order=0;end;
if nargin<7;OPT=1;end;
if nargin<6;Nmcs=100;end;
if nargin<5;method='autoscaling';end;
if nargin<4;K=10;end;
if nargin<3;A=2;end;


check=0; %+++ status variable:  1: Inf
[y,indexyy]=sort(y);
X=X(indexyy,:);
[Mx,Nx]=size(X);
Qs=floor(Mx*(1-1/K));

yytest=[];yp=[];
nLV=zeros(Nmcs,1);
predError=zeros(Nmcs,1);


for group=1:Nmcs
    perm=randperm(Mx);
    calk=perm(1:Qs);testk=perm(Qs+1:end);
    
    Xcal=X(calk,:);ycal=y(calk);
    Xtest=X(testk,:);ytest=y(testk);
    
    CV=plsldacv(Xcal,ycal,A,K,method,0,order);
    optPC=CV.optPC;
    if CV.check==1;check==1;break;end;
    LDA=plslda(Xcal,ycal,optPC,method);
    ypred=plsldaval(LDA,Xtest,ytest);
    yytest=[yytest;ytest];
    predError(group)=sum(sign(ypred)~=ytest)/length(ytest);
    yp=[yp;ypred;];
    nLV(group)=optPC;       
    if OPT==1;fprintf('The %d/%dth outer loop finished.\n',group,Nmcs);end;
end

%+++ Find the most frequently chosen nLV.
uniLV=unique(nLV);
for j=1:length(uniLV); freq(j)=length(find(nLV==uniLV(j)));end
[maxf,maxindex]=max(freq);
optPC=uniLV(maxindex(1));

%+++ output
if check==0
  F=roccurve(yytest,yp,0);
  error_rate=sum(sign(yp)~=yytest)/length(yp); 
  DCV.method=method;
  DCV.check=check;
  DCV.minCV=error_rate; 
  DCV.sensitivity=F.sensitivity;
  DCV.specificity=F.specificity;
  DCV.predError=predError;
  DCV.nLV=nLV;
  DCV.optPC=optPC;
elseif check==1
  DCV.method=method;
  DCV.check=check;  
end
  
  
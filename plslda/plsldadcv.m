function DCV=plsldadcv(X,y,A,K,method,OPT,order)
%+++ K-fold double cross validation Cross-validation for PLS-LDA
%+++ Input:  X: m x n  (Sample matrix)
%            y: m x 1  (measured property)
%            A: The max PC for cross-validation
%            K: fold. when K = m, it is leave-one-out CV
%       method: pretreatment method. Contains: autoscaling, center etc.
%          OPT: =1 Print process.
%               =0 No print.
%               pareto,minmax,center or none.
%+++ Order: =1  sorted,default. For CV partition.
%           =0  random. 
%+++ Output: Structural data: CV
%+++ Hongdong Li, Oct. 16, 2008.
%+++ Revised in Jan.12, 2009.

if nargin<7;order=1;end;
if nargin<6;OPT=1;end;
if nargin<5;method='autoscaling';end;
if nargin<4;K=10;end;
if nargin<3;A=2;end;


check=0; %+++ status variable:  1: Inf

if order==1
  [y,indexyy]=sort(y);
  X=X(indexyy,:);
else
  indexyy=randperm(length(y));
  X=X(indexyy,:);
  y=y(indexyy);
end


A=min([size(X,1)-ceil(length(y)/K) size(X,2) A]);
yytest=[];YR=[];
[Mx,Nx]=size(X);
groups = 1+rem(0:Mx-1,K);
yytest=[];yp=[];nLV=zeros(K,1);
for group=1:K
    testk = find(groups==group);  calk = find(groups~=group);
    Xcal=X(calk,:);ycal=y(calk);
    Xtest=X(testk,:);ytest=y(testk);
    
    CV=plsldacv(Xcal,ycal,A,K,method,0,order);
    if CV.check==1;check==1;break;end;
    LDA=plslda(Xcal,ycal,CV.optPC,method);
    ypred=plsldaval(LDA,Xtest,ytest);
    yytest=[yytest;ytest];
    yp=[yp;ypred;];
    nLV(group)=CV.optPC;       
    if OPT==1;fprintf('The %dth outer loop finished.\n',group);end;
end

%+++ Find the most frequently chosen nLV.
uniLV=unique(nLV);
for j=1:length(uniLV); freq(j)=length(find(nLV==uniLV(j)));end
[maxf,maxindex]=max(freq);
optPC=uniLV(maxindex(1));


%+++ output
if check==0
  F=roccurve(yytest,yp,0);
  error=sum(sign(yp)~=yytest)/Mx; 
  DCV.method=method;
  DCV.check=check;
  DCV.error=error; 
  DCV.Sensitivity=F.sensitivity;
  DCV.Specificity=F.specificity;
  DCV.nLV=nLV;
  DCV.optPC=optPC;
elseif check==1
  DCV.method=method;
  DCV.check=check;  
end
  
  
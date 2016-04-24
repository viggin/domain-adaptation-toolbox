function CV=plsldacv(X,y,A,K,method,OPT,order)
%+++ K-fold Cross-validation for PLS-LDA
%+++ Input:  X: m x n  (Sample matrix)
%            y: m x 1  (measured property)
%            A: The max PC for cross-validation
%            K: fold. when K = m, it is leave-one-out CV
%       method: pretreatment method. Contains: autoscaling, center etc.
%          OPT: =1 Print process.
%               =0 No print.
%               pareto,minmax,center or none.
%+++ Order: =1  sorted, default. For CV partition.
%           =0  random. 
%+++ Output: Structural data: CV
%+++ Hongdong Li, Oct. 16, 2008, lhdcsu@gmail.com
%+++ Revised in Jan.12, 2009.

if nargin<7,order=1;end
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
yytest=[];YR=zeros(Mx,A);
for group=1:K
    testk = find(groups==group);  calk = find(groups~=group);
    Xcal=X(calk,:);ycal=y(calk);
    Xtest=X(testk,:);ytest=y(testk);
    
    %data pretreatment    
    [Xcal,para1,para2]=pretreat(Xcal,method);
    if length(find(para2==0))>0; check=1;break;end %+++ Check data effectivness;
        
    ycals=pretreat(ycal,method);    
    Xtest=pretreat(Xtest,method,para1,para2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% SIM algorithm  %%%%%%%%%%%%%%
%     [B,C,P,T,U,wstar]=plssim(Xcal,ycals,A);    
    [B,wstar,T,P]=pls_nipals(Xcal,ycals,A,0);
    if sum(sum(isnan(T)))+sum(sum(isinf(T)))>0; check=1;break;end %+++ Check data effectivness;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:A
        %%%+++ train model.        
        TT=T(:,1:j); 
        C=ldapinv(TT,ycal,0);
%         coef=[wstar(:,1:j)*C(1:end-1);C(end)];        
        %+++ predict
%       y_est=Xtest*coef(1:end-1)+coef(end);
        Ttest=Xtest*wstar(:,1:j);
        y_est=Ttest*C(1:end-1)+C(end);
                
        YR(testk,j)=y_est;
    end
    if OPT==1;fprintf('The %dth fold for PLS-LDA finished.\n',group);end;
end

if check==0;

%+++ Original order
YR=YR(indexyy,:);
y=y(indexyy);
error=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:A
  error(i)=sum(sign(YR(:,i))~=y);    
end
error=error/Mx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mincv,index]=min(error);
ROC=roccurve(y,YR(:,index),0); %+++ ROC area.

k1=find(y==1);k2=find(y==-1);
y_est=sign(YR(:,index));
error_rate1=1-sum(y_est(k1)~=sign(y(k1)))/length(k1);
error_rate2=1-sum(y_est(k2)~=sign(y(k2)))/length(k2);

%+++ output  %%%%%%%%%%%%%%%%
CV.method=method;
CV.check=0;
CV.Ypred=YR;
CV.y=y;
CV.cv=error;
CV.minCV=mincv;
CV.optPC=index;
CV.sensitivity=error_rate1;
CV.specificity=error_rate2;
CV.AUC=ROC.AUC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else 
  CV.method=method;
  CV.check=1;      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+++ There is a song you like to sing.




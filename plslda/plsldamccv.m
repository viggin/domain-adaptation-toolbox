function CV=plsldamccv(X,y,A,method,N,ratio,OPT)
%+++ Monte Carlo Cross Validation for PLS-LDA
%+++ Input:  X: m x n  (Sample matrix)
%            y: m x 1  (measured property)
%            A: The max PC for cross-validation
%            N: The number of Monte Carlo Simulation.
%        ratio: The ratio of calibration samples to the total samples.
%       method: pretreatment method. Contains: autoscaling,
%               pareto,minmax,center or none.
%          OPT: =1 : print process.
%               =0 : don't print process.
%+++ Output: Structural data: CV
%+++ Hongdong Li, Oct. 16, 2008.
%+++ Revised in Jan.12, 2009.

if nargin<7;OPT=1;end;
if nargin<6;ratio=0.8;end
if nargin<5;N=1000;end
if nargin<4;method='autoscaling';end
if nargin<3;A=2;end



A=min([size(X,2) A]);
[Mx,Nx]=size(X);
Q=floor(Mx*ratio);
yytest=[];YR=[];
for i=1:N
    np=randperm(Mx);
    calk=np(1:Q);testk=np(Q+1:end);
    Xcal=X(calk,:);ycal=y(calk);
    Xtest=X(testk,:);ytest=y(testk);
    
    %data pretreatment    
    [Xcal,para1,para2]=pretreat(Xcal,method);
    ycals=pretreat(ycal,method);    
    Xtest=pretreat(Xtest,method,para1,para2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% SIM algorithm  %%%%%%%%%%%%%%
%     [B,C,P,T,U,wstar]=plssim(Xcal,ycals,A);    
    [B,wstar,T,P]=pls_nipals(Xcal,ycals,A,0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [B,W,T,P,Q,R2X,R2Y]=pls_nipals(X,Y,A,preproc)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yr=zeros(Mx-Q,A);
    for j=1:A
        %%%+++ train model.        
        TT=T(:,1:j);
        C=ldapinv(TT,ycal,0);
        coef=[wstar(:,1:j)*C(1:end-1);C(end)];        
        %+++ predict
        y_est=Xtest*coef(1:end-1)+coef(end);
        yr(:,j)=y_est;
    end
    YR=[YR;yr];yytest=[yytest;ytest];
    if OPT==1;fprintf('The %d/%dth sampling for MCCV PLS-LDA finished.\n',i,N);end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:A
  cv(i)=sum(sign(YR(:,i))~=yytest)/length(yytest);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mincv,index]=min(cv);
sespanalysis=sesp(yytest,YR(:,index));
%+++ MCCV test error
error_test=zeros(N,1);ntest=Mx-Q;
YRopt=YR(:,index);
for i=0:N-1
 ytestsub=yytest(i*ntest+1:i*ntest+ntest);
 ypresub=sign(YRopt(i*ntest+1:i*ntest+ntest));
 error_test(i+1)=sum(ytestsub~=ypresub)/ntest;    
end


% k1=find(y==1);k2=find(y==-1);
% y_est=sign(YR(:,index));
% error_rate1=1-sum(y_est(k1)~=sign(y(k1)))/length(k1);
% error_rate2=1-sum(y_est(k2)~=sign(y(k2)))/length(k2);

%+++ output  %%%%%%%%%%%%%%%%
CV.method=method;
CV.Ypred=YR;
CV.ytest=yytest;
CV.error_test=error_test;
CV.cv=cv;
CV.sensitivity=sespanalysis.sensitivity;
CV.specificity=sespanalysis.specificity;
CV.minCV=mincv;
CV.optPC=index;

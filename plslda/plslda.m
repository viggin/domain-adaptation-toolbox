function LDA=plslda(X,y,A,method)
%+++ programmed according to NIPALS algorithm by Hongdong Li,Oct. 2006.
%+++ A is the interation epoch.co is the regression coefficients linking x0 and y0
%+++ model: x=t*p'   y=t*r'=u*q'
%+++ y0 has to satisfy:  +1: positive class and -1:negative class;
%+++ method:    'none'
%        or 'center' 
%        or 'autoscaling';
%        or 'pareto';
%+++ Advisor£º Yizeng Liang, yizeng_liang@263.net
%+++ H.D. Li, Feb. 8, 2009, lhdcsu@gmail.com


if nargin<4;method='autoscaling';end;
if nargin<3;A=2;end;
if nargin<2;warn('Wrong input parameters!');end
A=min([size(X) A]);
check=0; %+++ check data effectiveness;

%%%%%%%%%% PLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xs,para1,para2]=pretreat(X,method);
[ys,ypara1,ypara2]=pretreat(y,method);    
[B,W,T,P,Q,R2X,R2Y]=pls_nipals(Xs,ys,A,0);
%%%%%%%%%% LDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if sum(sum(isnan(T)))+sum(sum(isinf(T)))>0; check=1;end %+++ Check data effectivness;

if check==0;
  [tpt,tpw,tpp,SR]=tp(Xs,B);             %+++ target projection and selectivity ratio
  C=ldapinv(T,y,0);                      %+++ Discriminant analysis
  yfit=T*C(1:end-1)+C(end);
  ROC=roccurve(y,yfit,0);
   %%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%
  LDA.method=method;
  LDA.scale_para=[para1;para2];
  LDA.yreal=y;
  LDA.regression_coef=B;
  LDA.weight=W;
  LDA.Xscores=T;
  LDA.Xloadings=P;
  LDA.R2X=R2X;
  LDA.R2Y=R2Y;
  LDA.tpScores=tpt;
  LDA.tpLoadings=tpp;
  LDA.SR=SR;
  LDA.LDA='*********  LDA ***********';
  LDA.check=0;
  LDA.coef_lda_pc=C;
  LDA.coef_lda_origin=[W(:,1:A)*C(1:end-1);C(end)];
  LDA.yfit=yfit;
  LDA.error=1-ROC.accuracy;
  LDA.sensitivity=ROC.sensitivity;
  LDA.specificity=ROC.specificity;
  LDA.AUC=ROC.AUC;
elseif check==1
  LDA.check=1; 
end

%+++ There is a song you like to sing.
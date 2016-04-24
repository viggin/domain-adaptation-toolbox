function F=scarsplslda(X,y,A,K,method,num,OPT) 
%+++ This is the simplified version of CARS, only retaining the EDF element.
%+++ num: the number of Monte Carlo Sampling.
%+++ A:   the maximal number of PLS components to extract.
%+++ fold: number of folds for cross validation.
%+++ method: pretreat method, 'center' or 'autoscaling'.
%+++ OPT: 1: Plot. 0: No plot.
%+++ Advisor: Yizeng Liang, yizeng_liang@263.net.
%+++ Hongdong Li, Jan.3, 2009.
%+++ Reference: Hongdong Li, Yizeng Liang, Qingsong Xu and Dongsheng Cao,
%         Key wavelengths screening using competitive adaptive reweighted sampling method
%         for multivariate calibration J?. Anal. Chim. Acta, 2009, 648 (1): 77-84 



%+++ Initial settings.
if nargin<7;OPT=1;end;  %+++ OPT==1: then figure output.
if nargin<6;num=50;end
if nargin<5;method='autoscaling';end
if nargin<4;K=5;end
if nargin<3;A=2;end

[Mx,Nx]=size(X);
A=min([Mx Nx A]);
index=1:Nx;


r0=1;
r1=2/Nx;
Vsel=1:Nx;

W=zeros(num,Nx);
Ratio=zeros(1,num);

%+++ Parameter of exponentially decreasing function. 
b=log(r0/r1)/(num-1);  a=r0*exp(b);

%+++ Main Loop
for iter=1:num     
     
     LDA=plslda(X(:,Vsel),y,A,method);    %+++ PLS model
     w=zeros(Nx,1);coef=LDA.coef_lda_origin(1:end-1);
     w(Vsel)=coef;W(iter,:)=w; 
     w=abs(w);                                  %+++ weights
     [ws,indexw]=sort(-w);                      %+++ sort weights
     
     ratio=a*exp(-b*(iter+1));                      %+++ Ratio of retained variables.
     Ratio(iter)=ratio;
     K=round(Nx*ratio);  
     
     
     w(indexw(K+1:end))=0;                      %+++ Eliminate some variables with small coefficients.  
     
     Vsel=find(w~=0);         
     fprintf('The %dth variable sampling finished.\n',iter);    %+++ Screen output.
 end

%+++  Cross-Validation to choose an optimal subset;
RMSEP=zeros(1,num);
Rpc=zeros(1,num);
for i=1:num
   vsel=find(W(i,:)~=0);
   CV=plsldacv(X(:,vsel),y,A,10,method,0);
   RMSEP(i)=min(CV.cv);
   Rpc(i)=CV.optPC;
   fprintf('The %dth subset finished.\n',i);
end
Rmin=min(RMSEP);
indexOPT=find(RMSEP==Rmin);
indexOPT=indexOPT(end);

%+++ output
F.method=method;
F.nPC=A;
F.W=W;
F.vim=sum(W,2);
F.cv=RMSEP;
F.minCV=Rmin;
F.iterOPT=indexOPT;
F.optPC=Rpc(indexOPT);
F.ratio=Ratio;
F.vsel=find(W(indexOPT,:)~=0)';
F.vim=abs(sum(W));

%+++ Plot
 if OPT==1;plotcars(F);end;
%+++ END



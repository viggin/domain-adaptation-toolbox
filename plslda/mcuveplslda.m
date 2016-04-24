function F=mcuveplslda(X,Y,A,K,method,N,ratio,Vmax,OPT)
%+++ uninformative variable elimination(UVE)-PLS.
%+++ Input:  X: m x n  (Sample matrix)
%            y: m x 1  (measured property)
%           PC: The max PC for cross-validation
%            N: The number of Monte Carlo Simulation.
%        ratio: The ratio of calibration samples to the total samples.
%       method: pretreatment method. Contains: autoscaling,
%               pareto,minmax,center or none.

%+++ Default parameters
if nargin<9;OPT=1;end
if nargin<8;Vmax=200;end;
if nargin<7;ratio=0.8;end
if nargin<6;N=1000;end
if nargin<5;method='autoscaling';end
if nargin<4;K=10;end
if nargin<3;A=2;end

%+++ Monte Carlo Uninformative Variable Elimination.
[Mx,Nx]=size(X);
Q=floor(Mx*ratio);
A=min([size(X) A]);
C=sparse(N,Nx);  

for group=1:N
      temp=randperm(Mx);
      calk=temp(1:Q);      
      Xcal=X(calk,:);ycal=Y(calk);  
      PLSLDA=plslda(Xcal,ycal,A,method);
      coef=PLSLDA.coef_lda_origin;
      C(group,:)=coef(1:end-1)';    
      if OPT==1; fprintf('The %dth sampling for MC-UVE-PLSLDA finished.\n',group);end
  end
  U=mean(C);  SD=std(C);  RI=abs(U./SD);
  [RIs,indexRI]=sort(-RI);
  Vsel=indexRI;
  
%+++ Variable evaluating
VR=[];
nLV=[];
for i=1:Vmax
    CV=plsldacv(X(:,Vsel(1:i)),Y,A,K,method,0);
    nLV(i)=CV.optPC;
    VR=[VR;[1-CV.minCV CV.sensitivity CV.specificity]];
    if OPT==1; fprintf('The %dth variable evaluating finished.\n',i); end
end

SD=100;Kopt=0;acc=0;
for i=1:size(VR,1)
   acctemp=VR(i,1);
   if acctemp>acc 
     acc=acctemp;  
     Kopt=i;
     SD=std(VR(i,:));
   elseif acctemp==acc
     tempSD=std(VR(i,:)); 
     if tempSD<SD; Kopt=i;SD=tempSD;end
   end
end

%+++ Output
F.RI=RI;
F.SortedVariable=Vsel;
F.Coefficient=C;
F.VariableEvaluation=VR;
F.nLV=nLV;
F.Kopt=Kopt;
F.BestVariables=Vsel(1:Kopt);
F.optPC=nLV(Kopt);
F.BestResults=VR(Kopt,:);


%+++ There is a song you like to sing.



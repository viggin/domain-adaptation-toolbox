%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  This script is used to test whether the functions in this  %%%%%%%%%%  
%%%%%%%%%%  package can run smoothly.                                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  H.D. Li, lhdcsu@gmail.com                                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%+++ Import data;
load DM2;
%+++ Cross validation
A=6;
K=5;
method='autoscaling';
N=500;
Nmcs=50;
CV=plsldacv(Xcal,ycal,A,K,method)
MCCV=plsldamccv(Xcal,ycal,A,method,N)
DCV=plsldardcv(Xcal,ycal,A,K,method,Nmcs)

%+++ Build a PLS-LDA model
nLV=3;
LDA=plslda(Xcal,ycal,nLV);
[ScoresTest]=plsldaproj(LDA,Xcal(1:3,:))
 
%+++ Scores plot
plotlda(LDA,1,0,[2 3 1]);
figure;
plotlda(LDA,1,1,[1  2 3]);
figure;
plotlda(LDA,0,0,[1 2]);
figure;
plotlda(LDA,1,1,[2 3 ]);
figure;
plotlda(LDA,2,1,[3 1]);

%+++ CARS-PLSLDA for variable selection
CARS=carsplslda(Xcal,ycal,A,K,method,50);
figure;
plotcars(CARS);
%+++ simplified version of CARS-PLSLDA for variable selection
sCARS=scarsplslda(Xcal,ycal,A,K,method,50);
figure;
plotcars(sCARS);

%+++ SPA for vairable selection: based on Model Population Analysis
SPA=spa(Xcal,ycal,A,K,method,N,0.7,15);
figure;
bar(SPA.COSS,'b','edgecolor','w');
figure;
plotspa(SPA,SPA.RankedVariable(1));
p=SPA.p(SPA.RankedVariable(1))
%+++ Test ended







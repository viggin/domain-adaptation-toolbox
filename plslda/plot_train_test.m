%+++ Import data;
clear;
load DM2;
np=randperm(96);
Xtrain=Xcal(np(1:60),:);
ytrain=ycal(np(1:60));
Xtest=Xcal(np(61:96),:);
ytest=ycal(np(61:96));


%+++ Cross validation
A=6;
K=5;
method='autoscaling';
%+++ Build a PLS-LDA model
nLV=3;
LDA=plslda(Xtrain,ytrain,nLV);
[ScoresTest]=plsldaproj(LDA,Xtest)

%+++ Plot
plotlda(LDA,0,1,[1 2 3],{'w','w','w'});            %+++ train set
% figure;
% plotlda(LDA,0,1,[1 2],{'w','w'})  %+++ test set
hold on;
kn=find(ytest==-1);
kp=find(ytest==1);
plot3(ScoresTest(kn,1),ScoresTest(kn,2),ScoresTest(kn,3),'rd');
plot3(ScoresTest(kp,1),ScoresTest(kp,2),ScoresTest(kp,3),'go');



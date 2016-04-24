function plotcars(F)
%%% F is the returned structural data obtained from arsplslda.m. 
%%% This function is to plot the result.
%+++ Hongdong Li, Jan 15, 2009.


W=F.W;
RMSECV=F.cv;
indexOPT=F.iterOPT;
num=length(RMSECV);

subplot(311);
for i=1:num;L(i)=length(find(W(i,:)~=0));end
plot(L,'linewidth',1.5);
xlabel('Number of MC samplings');ylabel('Number of sampled features');
subplot(312);
plot(1:num,RMSECV,'linewidth',1.5);
%   text(1:num,RMSEP,num2str(Rpc'));
xlabel('Number of MC samplings');ylabel('10-fold CVPA');

subplot(313);
plot(W);hold on;minW=min(min(W));maxW=max(max(W));plot(repmat(indexOPT,1,20),linspace(minW,maxW,20),'b*','linewidth',0.5);
d=abs(maxW-minW)*0.05;
axis([0 num minW-d maxW+d]);
xlabel('Number of MC samplings');ylabel('Model coefficients path');

%%%%%%%%+++ END +++%%%%%%%%


function [Fvalue,Vsel]=fishersel(X,y)
%+++ To perform the variable selection with Fisher principle of two classes
%+++ problem. The element of y has to be +1 or -1.
%+++ Apr.27 2008, Hongdong Li

%+++ rearrange the dataset according to y;
kp=find(y==1);kn=find(y~=1);
Xp=X(kp,:);Xn=X(kn,:);np=size(Xp,1);nn=size(Xn,1);
N=size(X,2);iter=1;Fvalue=zeros(1,N);
Meanw=mean(X);Meanp=mean(Xp);Meann=mean(Xn);
%+++ Fisher selecting
while iter<=N
  fenzi=(Meanw(iter)-Meanp(iter))^2+(Meanw(iter)-Meann(iter))^2;
  fenmu=sum((Xp(:,iter)-Meanp(iter)).^2)/(np-1)+ sum((Xn(:,iter)-Meann(iter)).^2)/(nn-1);
  Fvalue(iter)=fenzi/fenmu;    
  iter=iter+1;
end
Fvalue=Fvalue/max(Fvalue);
%+++ sort the Fvalue
[A,Vsel]=sort(-Fvalue);









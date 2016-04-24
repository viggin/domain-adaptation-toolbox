function C=ldapinv(X,y,flag);
%+++flag: =1    Bayes approximation.
%+++      =0    Fisher DA.
%+++ The last element in C is the bias term.
%+++Hongdong Li,Nov. 23.

A=length(y);B=length(find(y==1));C=A-B;
r1=A/B;r2=A/C;R=[];
kp=find(y==1);
kn=find(y==-1);
XX=[[X(kp,:) ones(B,1)];-[X(kn,:) ones(C,1)]];
R=[ones(B,1)*r1;ones(C,1)*r2];
BB=ones(A,1);
if flag==1;C=inv(XX'*XX)*XX'*BB;elseif flag==0; C=inv(XX'*XX)*XX'*R;end
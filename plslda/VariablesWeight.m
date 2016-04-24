function [Weight,Vsel]=VariablesWeight(X,y,thita)

k1=find(y==1);
k2=find(y~=1);
X1=X(k1,:);
X2=X(k2,:);

for i=1:size(X1,2)
[m1 n1]=size(X1(:,i));
[m2 n2]=size(X2(:,i));
m=m1+m2;
dis_12=dis(X1(:,i),X2(:,i));
dis_11=dis(X1(:,i),X1(:,i));
dis_22=dis(X2(:,i),X2(:,i));
W(i)=m1*m2/m^2*dis_12/(m1/m*dis_11+m2/m*dis_22);
% W(i)=(mean(X1(:,i))-mean(X2(:,i)))^2/(var(X1(:,i))+var(X2(:,i)));
end
Vsel=find(W>=thita);
Weight=W(Vsel);

function [Sdis]=dis(X1,X2)
[m1 n1]=size(X1);
[m2 n2]=size(X2);
XX1=repmat(X1,1,m2);
XX2=repmat(X2',m1,1);
XX=(XX1-XX2).^2;
Sdis=sum(sum(XX));

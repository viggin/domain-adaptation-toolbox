function t=tvalue(X,y)

k1=find(y==1);
k2=find(y~=1);

X1=X(k1,:);
X2=X(k2,:);
u1=mean(X1);
u2=mean(X2);
N1=length(k1);
N2=length(k2);
t=zeros(1,size(X,2));
for j=1:size(X,2);
    s1=sum((X1(:,j)-u1(j)).^2);
    s2=sum((X2(:,j)-u2(j)).^2);
    sigma =sqrt((s1+s2)/(N1+N2-2));
    sej=sigma/sqrt(1/N1+1/N2);
    t(j)=(u2(j)-u1(j))/sej;
end









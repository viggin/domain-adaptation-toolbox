function F=tstatis(X,y,P,N)

if nargin<4;N=3;end
if nargin<3;P=1000;end

[M Nvar]=size(X);
t0=tstatis1(X,y);
k1=find(y==1);
k2=find(y~=1);
T=zeros(P,Nvar);
TMC=zeros(N,Nvar);
Q=floor(M*0.9);


%++++ Reliability index
i=1;
while i<=N
  index=randperm(M);   % MCS
  cal=index(1:Q);
  tmc=tvalue(X(cal,:),y(cal));
  TMC(i,:)=tmc;
  fprintf('The %dth MCS finished.\n',i);
  i=i+1;
end   
RI=abs(mean(TMC)./std(TMC));



%+++ Permutation test for p value computation.
i=1;
while i<=P
  index=randperm(M);   % permutation
  yr=y(index);
  if sum(abs(yr-y))==0;continue;end
  tr=tvalue(X,yr);
  T(i,:)=tr;
  fprintf('The %dth permutation finished.\n',i);
  i=i+1;
end    

p=zeros(1,Nvar);
for i=1:Nvar
   p(i)=length(find(abs(T(:,i))>abs(t0(i))))/P;    
end


%+++ Output
F.t0=abs(t0)/max(abs(t0));
F.TMC=TMC;
F.RI=RI;
F.tr=T;
F.p=p;


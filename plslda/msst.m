function F=msst(X,y,N)
%+++ Multi-Scale significance test
if nargin<3;N=1000;end


[n,p]=size(X);
r0=15/p;
ratio=[0.5 0.6 0.7 0.8 0.9 0.95];
for i=1:N
  for j=1:length(ratio)
     
    [Xcal,ycal]=traintestselect(X,y,ratio(j));
    kp=find(ycal==1);
    kn=find(ycal~=1);
    for k=1:p
%      [ptemp,h] = ranksum(Xcal(kp,k),Xcal(kn,k));
     [h,ptemp] = ttest2(Xcal(kp,k),Xcal(kn,k));
     P(i,j,k)=-log10(ptemp);
    end
    
  end
  fprintf('The %ith sampling finished.\n',i);
end
%+++ Output
F.P=P;

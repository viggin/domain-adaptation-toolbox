function MTC=mtc(p,alpha,w)
%+++ Multiple Testing Correction for controling familywise error rate(FWER)using Holm-Bonferroni method.
%+++ w: weight, defaul, equal weight.
%+++ alpha: significance level of FWER;
%+++ p: unajusted p-value resulting from each single comparison/test.
%+++ Advisor£º Yizeng Liang, yizeng_liang@263.net
%+++ H.D. Li, Apr. 18, 2010, lhdcsu@gmail.com


if nargin<3;N=length(p);w=repmat(1/N,1,N);end
if nargin<2;alpha=0.05;end


%+++ Sorting p-value from low to high.
[sp,index]=sort(p);
sw=w(index);
reject=zeros(1,N);
ajustedp=zeros(1,N);


%+++ whether reject and computing ajusted p value
BH=zeros(N,2);
ajustedp(1)=min([1, sp(1)/sw(1)]);
flag=0;
for i=1:N
  spi=sp(i);
  threshold=alpha/(N-i+1);
  if i>1;ajustedp(i)=min([1, max([ajustedp(i-1) sum(sw(i:end))*sp(i)/sw(i)])]);end  
  
  if spi<threshold;
      if flag==0
        reject(i)=1;      
      end 
  else
      flag=1;
  end
  BH(i,:)=[spi threshold];
end

BH(index,:)=BH;
ajustedp(index)=ajustedp;
reject(index)=reject;
%+++ Output
MTC.BH=BH;
MTC.ajustedp=ajustedp;
MTC.reject=reject;

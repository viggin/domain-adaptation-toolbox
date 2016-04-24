function sel=weightsampling(w)
%+++ Bootstrap sampling
%+++ w: sampling weight, must be positive numbers.
%+++ Ns: Number of objects(samples or variables or others) to be sampled.
%+++ 2007.9.6,H.D. Li.

w=w/sum(w);
N1=length(w);
min_sec(1)=0; max_sec(1)=w(1);
for j=2:N1
   max_sec(j)=sum(w(1:j));
   min_sec(j)=sum(w(1:j-1));
end
% figure;plot(max_sec,'r');hold on;plot(min_sec);
sel=zeros(1,N1); 
for i=1:N1
  bb=rand(1);
  ii=1;
  while (min_sec(ii)>=bb | bb>max_sec(ii)) & ii<N1;
    ii=ii+1;
  end
  sel(i)=ii;
end      % w is related to the bootstrap chance







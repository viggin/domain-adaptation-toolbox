function vip=vipp(x,y,t,w)
%+++ vip=vipp(x,y,t,w);
%+++ t: socres, which can be obtained by pls_nipals.m 
%+++ w: weight, which can be obtained by pls_nipals.m 
%+++ to calculate the vip for each variable to the response;
%+++ vip=sqrt(p*q/s);

%initializing
[m,p]=size(x);
[m,h]=size(t);
[p,h]=size(w);
%calculate s
for i=1:h
    corr=corrcoef(y,t(:,i));
    co(i,1)=corr(1,2)^2;
end
s=sum(co);
%calculate q;
for i=1:p
    for j=1:h
        d(j,1)=co(j,1)*w(i,j)^2;
    end
    q=sum(d);
    vip(i,1)=sqrt(p*q/s);
end


    

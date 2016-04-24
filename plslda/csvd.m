function  [U,S,V]=csvd(X)
%+++ csvd means clever svd. Hehe...
%+++ Hongdong Li,Jul.28,2008
   
[Mx,Nx]=size(X);
if Mx<=Nx
    [u1,s1,v1]=svd(X*X',0);
     U=u1;S=sqrt(s1);V=X'*U*(inv(S))';
else
    [u1,s1,v1]=svd(X'*X,0);
    V=v1;S=sqrt(s1);U=X*V*inv(S);
end




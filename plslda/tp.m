function [t,w,p,sr]=tp(X,b)
%+++ target projection for PLS
%+++ regression vector from PLS


w=b/norm(b);
t=X*w;
p=X'*t/(t'*t);

Xtp=t*p';
Xr=X-Xtp;

%+++ Compute selectivity ratio
for i=1:size(X,2)
  vart(i)=sumsqr(X(:,i));
  vartp(i)=sumsqr(Xtp(:,i));
  varr(i)=sumsqr(Xr(:,i));
end
sr=vartp./(varr+eps); % sr.
%+++ end
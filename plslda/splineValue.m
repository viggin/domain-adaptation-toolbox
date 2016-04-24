function [xnew,ynew]=splineValue(x,y)
%+++ spline for interpolation.
%+++ 
%+++ Aug., 2010.

minx=min(x);
maxx=max(x);
d=(maxx-minx)/128;
xnew=linspace(minx-d,maxx+d,128)';
ynew=spline(x,y,xnew);



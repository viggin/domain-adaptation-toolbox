function classplot2random(X,y,flag,colorstr,MarkerSize)
%+++ flag: 0 No marker;
%          1 using clorstr
%          2 using class number;      

if nargin<5;MarkerSize=5;end
if nargin<4;colorstr={'g.';'b.';'k*';'gp';'cs';'g>';'c.';};end
if nargin<3;flag=0;end

nclass=unique(y);
n=length(y);
nperm=randperm(n);
hold on;
for i=1:n
  permi=nperm(i);  
  if y(permi)==nclass(1);colortemp=colorstr{1};else; colortemp=colorstr{2};end
  plot(X(permi,1),X(permi,2),colortemp,'MarkerSize',MarkerSize);    
end


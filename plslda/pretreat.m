function [X,para1,para2,Problem]=pretreat(X,method,para1,para2)
%   data pretreatment
Problem=0;
if nargin==2
  [Mx,Nx]=size(X);
   if strcmp(method,'autoscaling')
    para1=mean(X);para2=std(X);
   elseif strcmp(method,'center')
    para1=mean(X);para2=ones(1,Nx);
   elseif strcmp(method,'minmax')
    para1=min(X);maxv=max(X);para2=maxv-para1;  
   elseif strcmp(method,'pareto')
    para1=mean(X);para2=sqrt(std(X));
   elseif strcmp(method,'none')
    para1=zeros(1,Nx);para2=ones(1,Nx);
   else
    display('Wrong data pretreat method!');
   end
   
   if length(find(para2==0))<1
     for i=1:Nx; X(:,i)=(X(:,i)-para1(i))/para2(i);end
   else
     Problem=1;  
   end
   
   
   
elseif nargin==4
   [Mx,Nx]=size(X);
   for i=1:Nx     
     X(:,i)=(X(:,i)-para1(i))/para2(i);
   end
end

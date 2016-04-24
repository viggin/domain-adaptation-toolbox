function [X,index]=rmvconstvar(X,threshold)


[n,p]=size(X);
index=[];
for i=1:p
    L=length(unique(X(:,i)));
    r=1-L/n;
    if r>=threshold;index=[index;i];end
end
X(:,index)=[];


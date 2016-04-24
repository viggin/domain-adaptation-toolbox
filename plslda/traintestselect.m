function [Xcal,ycal,Xtest,ytest]=traintestselect(X,y,ratio,OPT)
%+++ This function is used to select a training and corresponding test
%    from the whole dataset.
%+++ X: sample matrix of size m x p
%+++ y: class label vector of size m x 1, must of element 1 or -1
%+++ ratio: The ratio of the number of train samples to the whole data
%+++ OPT: Selection methods
%      0: randomly selection
%      1: Kennard-Stone selection
%+++ Advisor: Yi-Zeng Liang, yizeng_liang@263.net
%+++ H.D. Li, Apr. 21, 2010, lhdcsu@gmail.com


%+++ Initial settings
if nargin<4;OPT=0;end
if nargin<3;ratio=0.8;end
if nargin<2;warning('Not enough input parameters!\n');end


%+++ ID number of two classes of samples
id1=find(y==1);
id2=find(y==-1);
n1=length(id1);
n2=length(id2);
X1=X(id1,:);X2=X(id2,:);
y1=y(id1,:);y2=y(id2,:);

%+++ Selection strategy
if OPT==0
  rank1=randperm(n1);
  rank2=randperm(n2);
elseif OPT==1
  rank1=ks(X1);
  rank2=ks(X2);
end

%+++ Number of calibration samples for each class
nc1=round(n1*ratio);
nc2=round(n2*ratio);

calk1=rank1(1:nc1);
calk2=rank2(1:nc2);
testk1=rank1((1+nc1):n1);
testk2=rank2((1+nc2):n2);

%+++ Output
Xcal=[X1(calk1,:);X2(calk2,:)];
ycal=[y1(calk1,:);y2(calk2,:)];

Xtest=[X1(testk1,:);X2(testk2,:)];
ytest=[y1(testk1,:);y2(testk2,:)];


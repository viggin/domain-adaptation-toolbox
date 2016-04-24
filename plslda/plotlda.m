function plotlda(plsmodel,flag,LINE,DimIndex,colorstr)
%+++ plsmodel: PLS-LDA model obtained by plslda.m
%+++ flag: Plot types. flag={0,1,2,3,4}; 
%+++ LINE: whether showing the discriminant line, default: no
%+++ DimIndex: The latent variables to be ploted,default DimIndex=[1 2];

if nargin<5; colorstr={'rd';'bo';'k*';'gp';'cs';'g>';'c.';};end
if nargin<4;DimIndex=[1 2];end
if nargin<3;LINE=0;end
if nargin<2;flag=0;end


y=plsmodel.yreal;
t=plsmodel.Xscores;

if max(DimIndex)>size(t,2);DimIndex=[1 2];end

Dim=length(DimIndex);

%+++ Plot
if Dim==2
  t11=t(:,DimIndex(1));t22=t(:,DimIndex(2));
  c=ldapinv(t(:,DimIndex),y,0); 
  t1=linspace(min(t11),max(t11),20);
  t2=-(t1*c(1)+c(3))/c(2);
  classplot2(t(:,DimIndex),y,flag,colorstr);hold on;
  if LINE==1;plot(t1,t2,'k-');end;
  d1=(max(t1)-min(t1))/10; d2=(max(t22)-min(t22))/10;
  labelx=sprintf('PLS-%d',DimIndex(1));
  labely=sprintf('PLS-%d',DimIndex(2));  
  xlabel(labelx);ylabel(labely);
  
elseif Dim==3
  len=25;
%   c=plsmodel.coef_lda_pc;
  c=ldapinv(t(:,DimIndex),y,0);
  
  t1=t(:,DimIndex(1));t2=t(:,DimIndex(2));t3=t(:,DimIndex(3));
  a1=min(t1);b1=max(t1);c1=min(t2);d1=max(t2);
  t11=linspace(a1,b1,len)';t22=linspace(c1,d1,len)';
  [T1,T2]=meshgrid(t11,t22);
  T3=-(T1*c(1)+T2*c(2)+c(4))/c(3);
  colormap([1  1  0; 0  1  1]);
  
  if LINE==1;surf(T1,T2,T3,ones(len,len)*0.3);hold on;end;  
  classplot2(t(:,DimIndex),y,flag,colorstr); 
  
  shading flat;
  d1=(max(t1)-min(t1))/50; d2=(max(t2)-min(t2))/50;d3=(max(max(T3))-min(min(T3)))/50;
  axis([min(t1)-d1 max(t1)+d1 min(t2)-d2 max(t2)+d2 min(min(T3))-d3 max(max(T3))+d3]);
  
  labelx=sprintf('PLS-%d',DimIndex(1));
  labely=sprintf('PLS-%d',DimIndex(2));  
  labelz=sprintf('PLS-%d',DimIndex(3));  
  xlabel(labelx);ylabel(labely);zlabel(labelz);
 
  
else
  disp('Dimension should be 2 or 3!');
end

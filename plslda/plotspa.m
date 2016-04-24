function [Q2in,Q2out]=plotspa(F,index,n0,n1)
%+++ plot the predictive accuracy distribution for a given variable
%+++ advisor: Yi-Zeng Liang, yizeng_liang@263.net
%+++ Coder: Hong-Dong Li, lhdcsu@gmail.com
%+++ Sep. 12, 2009

if nargin<4;n1=25;end
if nargin<3;n0=20;end;
if nargin<2;index=1;end


c=[1 1 1];  %+++ white color;
error0=F.error0;
error1=F.error1;
k=find(~isnan(error1(:,index))==1);
Q2in=error0(k);
errori=error1(:,index);
Q2out=errori(k);
% [n1,xout1] = hist(Q2in);
% [n0,xout0] = hist(Q2out);
% % axiss=[min([xout0 xout1])  max([xout0 xout1]) min([n0 n1])   max([n0 n1])];

c0=[1 1 1];
c1=[0.5 0.5 0.5];
[H,xbin1,ybin1]=histfitnew(Q2in(:,1),n1);
h = get(H,'Children');set(h,'FaceColor','g');
hold on;
[H,xbin0,ybin0]=histfitnew(Q2out(:,1),n0);box on;
h = get(H,'Children');set(h,'FaceColor','b');
xlabel('Prediction error');
ylabel('Probability density');


%+++ Axis limitation
minx=min([xbin0 xbin1]);
maxx=max([xbin0 xbin1]);
dx=(maxx-minx)/20;
miny=min([ybin0 ybin1]);
maxy=max([ybin0 ybin1]);
dy=(maxy-miny)/20;
axis([minx-dx maxx+dx miny-0.01*dy maxy+dy]);
box on;



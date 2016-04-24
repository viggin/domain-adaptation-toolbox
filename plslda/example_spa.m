%+++ Data load
load T2DM
%+++ Subwindow Permutation Analysis
F=spa(X,y,5,5,'autoscaling',300,0.8,10);
bar(F.COSS);
xlabel('Variable index');
figure;
plotspa(F,11);
vargrouping(F,[1:21],1)




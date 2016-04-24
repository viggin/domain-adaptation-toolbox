function plot2D(X,Y)
% {plot2D} plots a binary classification bidimensional dataset.
%
%      plot2D(X,Y)
%      
%      X: N-by-2 matrix of N examples 
%      Y: vector with N targets in {-1,0,+1}, where 0 means 'unlabeled'
% 
%      D: M-by-N distance matrix
%
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it
%         * based on the code of Vikas Sindhwani, vikas.sindhwani@gmail.com

pos=(Y==1);
neg=(Y==-1);
unlab=(Y==0);

rmin=min(min(X))-0.2;
rmax=max(max(X))+0.2;
xlim([rmin,rmax]);
ylim([rmin,rmax]);

plot(X(unlab,1),X(unlab,2),'ks','MarkerSize',3);
hold on;  
plot(X(pos,1),X(pos,2),'rd','MarkerSize',8, ...
     'MarkerFaceColor','r','MarkerEdgeColor','k'); 
hold on;
plot(X(neg,1),X(neg,2),'bo' ,'MarkerSize',8, ...
    'MarkerFaceColor','b','MarkerEdgeColor','k');
hold on;


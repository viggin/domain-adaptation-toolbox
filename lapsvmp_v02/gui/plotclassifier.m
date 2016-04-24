function plotclassifier(classifier,xrange,yrange)
% {plotclassifier} plots the output of a trained classfier on a 2D grid.
%
%      plotclassifier(classifier,xrange,yrange)
%      
%      classifier: the trained classifier structure (generated with the
%                  'saveclassifier' function)
%      xrange: [xmin, xmax]
%      yrange: [ymin, ymax]
%
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it
%         * based on the code of Vikas Sindhwani, vikas.sindhwani@gmail.com

% getting options
x=xrange;
y=yrange;
alpha=classifier.alpha;
xtrain=classifier.xtrain;
b=classifier.b;
name=classifier.name;
kerneltype=classifier.options.Kernel;
kernelparam=classifier.options.KernelParam;
lambdas=[classifier.options.gamma_A, classifier.options.gamma_I];

% making a grid of points
[xx,yy]=meshgrid(x,y);
X=[xx(:) yy(:)];

% computing outputs on the grid
K=calckernel(classifier.options, xtrain, X);
z=real(K*alpha + b);
z(z>1)=1;
z(z<-1)=-1;

% drawing
Z=reshape(z,length(x),length(y));
contourf(xx,yy,Z,[0 0 0],'LineStyle','none');

hold on;
title([name '  ('  kerneltype  ', \sigma =' num2str(kernelparam) ')']);
if ~isempty(strfind(name,'lap'))
    xlabel([' \gamma_A = ' num2str(lambdas(1)) ...
            ' \gamma_I = ' num2str(lambdas(2))]);
else
    xlabel([' \gamma_A = ' num2str(lambdas(1))]);
end

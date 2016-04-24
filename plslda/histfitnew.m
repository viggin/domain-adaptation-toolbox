function [H,xbin,ybin]= histfitnew(data,nbins)
%HISTFIT Histogram with superimposed fitted normal density.
%   HISTFIT(DATA,NBINS) plots a histogram of the values in the vector DATA.
%   using NBINS bars in the histogram. With one input argument, NBINS is set 
%   to the square root of the number of elements in DATA. 
%
%   H = HISTFIT(DATA,NBINS) returns a vector of handles to the plotted lines.
%   H(1) is a handle to the histogram, H(2) is a handle to the density curve.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 2.13.2.1 $  $Date: 2003/11/01 04:26:28 $
%

[row,col] = size(data);
if min(row,col) > 1
   error('stats:histfit:VectorRequired','DATA must be a vector.');
end

if row == 1
  data = data(:);
end
row = sum(~isnan(data));

if nargin < 2
  nbins = ceil(sqrt(row));
end
band=(max(data)-min(data))/nbins;
[n,xbin]=hist(data,nbins);

mr = nanmean(data); % Estimates the parameter, MU, of the normal distribution.
sr = nanstd(data);  % Estimates the parameter, SIGMA, of the normal distribution.

x=(-3*sr+mr:0.1*sr:3*sr+mr)';% Evenly spaced samples of the expected data range.

% hh = bar(xbin,n/sum(n)/band,1); % Plots the histogram.  No gap between bars.
% np = get(gca,'NextPlot');    
% set(gca,'NextPlot','add')    
                             
xd = data;

rangex = max(xd(:)) - min(xd(:)); % Finds the range of this data.
binwidth = rangex/nbins;    % Finds the width of each bin.


y = normpdf(x,mr,sr);  
% y = row*(y*binwidth);   % Normalization necessary to overplot the histogram.

% hh1 = plot(x,y,'r-','LineWidth',2);     % Plots density line over histogram.
hold on;
ybin=max(y)*n/max(n);
H= bar(xbin,ybin,1,'EdgeColor','w');           % Plots the histogram.  No gap between bars.
hh1 = plot(x,y,'r-','LineWidth',1.5);     % Plots density line over histogram.
% if nargout == 1
%   h = [hh; hh1];
% end
% 
% set(gca,'NextPlot',np) 

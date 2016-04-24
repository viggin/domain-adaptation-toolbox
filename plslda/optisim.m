function [subset]=optiSim(x,k,S,fig)

%
% Function: [subset]=OptiSim(x,k,S,fig)
%
% Aim: uniform selection of objects
%
% Input: x- data
%        k- number of objects in the subset to be selected
%        S- number of objects in the subsample 
%        fig- plot 1/Yes, 0/No (feature 1 vs. feature 2)
%
% Output: subset- indices of objects selected from x to subset
%
% Date: 15.10.01
% Author: Michal Daszykowski
%
% Copyright(c) 2002 for ChemoAC
% FABI, Vrije Universiteit Brussels
% Laarbeeklaan 103, 1090 Brussels	 
%
% References: R. D. Clark, OptiSim: An extended Dissimilarity 
% Selection method for finding diverse representative subsets,
% J. Chem. Inf. Comput. Sci. 37 (1998) 118-1188

[m,n]=size(x);
a=randperm(m);												% mix the objects in x

x=[[1:m]' x(a,:)];											% add indices

xx=x;														% rename oryginal data
[m,n]=size(x);

D=dist(mean(x(:,2:n)),x(:,2:n));							% calculate the distance to the mean
[i1 ind]=min(D);			
object=x(ind,:);											% frind the closest object to the mean
x(ind,:)=[];												% remove this object from x

bin=[];														% set recycle bin as empty
subsample=[];												% set subsample as empty
subset=object(1);											% put the object into subset

% Thresold													% calculate the treshold

Eps=((prod(max(x(:,2:n))-min(x(:,2:n)))*m/k*gamma(.5*(n-1)+1))/(m*sqrt(pi.^(n-1)))).^(1/(n-1));

while size(subset,1)<k   
      
	while size(subsample,1)<S

        if isempty(x);break;end								% stop if there are no objects in x

		object=x(1,:);										% select the 1st object
		D=dist(object(2:n),xx(subset,2:n));					% calculate the distance to objects in subset
             
		if min(D)<Eps 										% if the distance is lower than Eps
			x(1,:)=[];										% remove the object from x
			D=[];
		else 
			subsample=[subsample;[object(1) min(D)]];		% add the object to the subsample 
			x(1,:)=[];								 		% remove an object from x
		end     	                       
	end 													% <-- end of the 2nd while 

	if isempty(x);
		x=[x;xx(bin,:)];									% if x is empty add the bin to x 
		bin=[];
	end     

	if size(subsample,1)<1;break;end						% stop if subsample is empty	
	
	[i1 i2]=max(subsample(:,2));							% find an object with the highest D
 	subset=[subset;subsample(i2,1)];						% add it to the subsample
	subsample(i2,:)=[];										% remove the object from subsample
	bin=[bin;subsample(:,1)];								% put the remaining objects to the bin
	subsample=[]; 											% clear subsample
end															% <-- end of the 1st while


if fig==1													% plot figure
	plot(xx(:,2),xx(:,3),'.','markeredgecolor',[.4 .4 .4])
	hold on
	plot(xx(subset,2),xx(subset,3),'ks','markerfacecolor','k')
	xlabel('feature 1')
	ylabel('feature 2')
	figure(gcf);
end

subset=a(subset);											% sort the subset

 
%....................................................
function [D]=dist(i,x)

% Aim: calculates Euclidean distances between i and x	 

D=sqrt(sum((((ones(size(x,1),1)*i)-x).^2)'));


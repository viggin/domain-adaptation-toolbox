function classplot2(X,y,flag,colorstr)
%+++ flag:  0 No marker;
%           1 using clorstr
%           2 using class number;     
%+++ Advisor: YZ Liang, yizeng_liang@263.net
%+++ H.D. Li, lhdcsu@gmail.com

if nargin<4; colorstr={'rd';'bo';'k*';'gp';'cs';'g>';'c.';};end
nclass=unique(y);
if nclass>length(colorstr); warning('Not enough items for the 4th input parameters!');end

d=(max(X(:,1))-min(X(:,1)))/50;
if size(X,2)==2
 iter=1;
 for i=1:length(nclass)
   index_i=find(y==nclass(i));
   Xsub=X(index_i,:);
   if flag==0
     plot(Xsub(:,1),Xsub(:,2),colorstr{i});   
   elseif flag==1;
     plot(Xsub(:,1),Xsub(:,2),colorstr{i});
     text(Xsub(:,1)+d,Xsub(:,2),num2str(iter));
   elseif flag==2;    
     plot(Xsub(:,1),Xsub(:,2),colorstr{i},'marker','none');  
     text(Xsub(:,1)+d,Xsub(:,2),num2str(iter));  
   elseif flag==3     
     plot(Xsub(:,1),Xsub(:,2),colorstr{i});
     text(Xsub(:,1)+d,Xsub(:,2),num2str([index_i]));    
   elseif flag==4     
     plot(Xsub(:,1),Xsub(:,2),colorstr{i});
     text(Xsub(:,1)+d,Xsub(:,2),num2str([1:size(Xsub,1)]')); 
   end
   hold on;   
   iter=iter+1;
 end
elseif size(X,2)==3
 iter=1;
 for i=1:length(nclass)
   index_i=find(y==nclass(i));
   Xsub=X(index_i,:);
   if flag==0
     plot3(Xsub(:,1),Xsub(:,2),Xsub(:,3),colorstr{i});   
   elseif flag==1;
     plot3(Xsub(:,1),Xsub(:,2),Xsub(:,3),colorstr{i});
     text(Xsub(:,1),Xsub(:,2),Xsub(:,3),num2str(iter));
   elseif flag==2;    
     plot3(Xsub(:,1),Xsub(:,2),Xsub(:,3),colorstr{i},'marker','none');  
     text(Xsub(:,1),Xsub(:,2),Xsub(:,3),num2str(iter));  
   elseif flag==3     
     plot3(Xsub(:,1),Xsub(:,2),Xsub(:,3),colorstr{i});
     text(Xsub(:,1),Xsub(:,2),Xsub(:,3),num2str([index_i]));    
   elseif flag==4     
     plot3(Xsub(:,1),Xsub(:,2),Xsub(:,3),colorstr{i});
     text(Xsub(:,1),Xsub(:,2),Xsub(:,3),num2str([1:size(Xsub,1)]')); 
   end
   hold on;   
   iter=iter+1;
 end    
else
  warning('Wrong input parameters! Dimension is larger than 3.');    
end








    


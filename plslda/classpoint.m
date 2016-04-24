function [X1,X2]=classpoint
X1=[];
X2=[];
axis([0 1 0 1]);hold on;
button=1;
while button==1
[x,y,button] = ginput(1);
plot(x,y,'b+');
X1=[X1;[x y]];
end
button=1;
while button==1
[x,y,button] = ginput(1);
plot(x,y,'r*');
X2=[X2;[x y]];
end
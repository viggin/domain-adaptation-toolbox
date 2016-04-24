function result=sesp(yreal,ypred)
%+++ To calculate the sensitivity(Se) and the specificity(Sp) of 
%+++ a binary classification problem.
%+++ y-value has to be 1 or -1.
%+++ Hongdong Li, Apr.29,2008.

ypred=sign(ypred);
yreal=sign(yreal);
p=0;n=0;o=0;u=0;LEN=length(yreal);
for i=1:LEN
    if  yreal(i)==1
       if ypred(i)==1;p=p+1;else;o=o+1;end
    elseif yreal(i)~=1
       if ypred(i)~=1; n=n+1;else;u=u+1;end
    end   
end

%+++ output
result.sensitivity=p/(p+o);
result.specificity=n/(n+u);
result.accuracy=(p+n)/LEN;




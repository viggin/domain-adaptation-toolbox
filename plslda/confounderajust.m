function F=confounderajust(XC,XI,confounder,OPT)
%+++ AIM: Remove the variance caused by confunding variables by using least squares.
%+++ XC: Independent variable matrix including confounding factors.
%+++ XI: The interesting variable data to be corrected.
%+++ confounder: the column index of XC chosen as confunding factors.


[Mcal,Ncal]=size(XC);

%+++ Partition data and linear modelling between XI and XC
Interest=1:Ncal;
Interest(confounder)=[];
XCi=[XC(:,Interest) ones(Mcal,1)];
XCp=XC(:,confounder);
Bin=inv(XCi'*XCi)*XCi'*XCp;
E=XCp-XCi*Bin;
R2b=1-sumsqr(E)/sumsqr(XCp);

%+++ Linear modeling correlating y and X
XC=[XC ones(Mcal,1)];
B=inv(XC'*XC)*XC'*XI;

Bc=B(confounder,:);
if OPT==0
  XIc=XI-XC(:,confounder)*Bc;
elseif OPT==1
  XIc=XI-E*Bc;
end
MEAN=mean(XI);
R2=sumsqr(XIc)/sumsqr(XI);
%+++ Output
F.confounder=confounder;
F.B=B;
F.Bin=Bin;
F.R2=R2;
F.R2b=R2b;
F.XI_ajusted=XIc;










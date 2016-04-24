%# Function [B,Wstar,T,P,Q,R2X,R2Y,W]=pls_nipals(X,Y,A,preproc);
%#
%# AIM:         performs PLS calibration on X and Y
%# PRINCIPLE:   Uses the NIPALS algorithm to perform PLS model calibration
%# REFERENCE:   Multivariate Calibration, H. Martens, T. Naes, Wiley and
%#              sons, 1989
%#
%# INPUT:
%# X            matrix of independent variables (e.g. spectra) (n x p)
%# Y            vector of y reference values (n x 1)
%# A            number of PLS factors to consider
%# preproc      preprocessing applied to data
%#              0: no preprocessing
%#              1: column mean-centering of X and Y
%#
%# OUTPUT:
%# B            regression coefficients (p x 1)
%# W            X-weights (p x A)
%# T            scores (n x A)
%# P            X-loadings (p x A)
%# Q            Y-loadings (A x 1)
%# R2X          percentage of X variance explained by each PLS factor
%# R2Y          percentage of Y-variance explained by each PLS factor
%#
%# AUTHOR:      Xavier Capron
%# 			    Copyright(c) 2004 for ChemoAC
%# 			    FABI, Vrije Universiteit Brussel
%# 			    Laarbeeklaan 103, 1090 Jette
%# 			    Belgium
%#             
%# VERSION: 1.0 (24/11/2004)


function [B,Wstar,T,P,Q,R2X,R2Y,W]=pls_nipals(X,Y,A,preproc)

[n,p]=size(X);

if preproc==1
    [X,mX]=center(X);
    [Y,mY]=center(Y);
end

Xorig=X;
Yorig=Y;

ssqX=sum(sum((X.^2)));
ssqY=sum(Y.^2);

for a=1:A    
    W(:,a)=X'*Y;
    W(:,a)=W(:,a)/norm(W(:,a));
    T(:,a)=X*W(:,a);
    P(:,a)=X'*T(:,a)/(T(:,a)'*T(:,a));
    Q(a,1)=Y'*T(:,a)/(T(:,a)'*T(:,a));
    X=X-T(:,a)*P(:,a)';
    Y=Y-T(:,a)*Q(a,1);
    R2X(a,1)=(T(:,a)'*T(:,a))*(P(:,a)'*P(:,a))/ssqX*100;
    R2Y(a,1)=(T(:,a)'*T(:,a))*(Q(a,1)'*Q(a,1))/ssqY*100;   
end
Wstar=W*(P'*W)^(-1);
B=Wstar*Q;




%# function [B,C,P,T,U,R,R2X,R2Y]=plssim(X,Y,A,S,XtX)			    
%#									    
%# AIM:	      Full implementation of SIMPLS approach to PLS regression	     
%#	      for (multivariate) Y.					    
%#									    
%# PRINCIPLE: Ref.: S.de Jong, Chemom.Intell.Lab.Syst.,18 (1993) 251-263.   
%#									    
%#									    
%# INPUT:     X    n  x px   predictor  data               <centered>	    
%#            Y    n  x m    predictand data               <centered>	    
%#	      A    1  x 1    max # factors to consider     <optional>	    
%#	      S    px x m    X'*Y                          <optional>	    
%#	      XtX  px x px   X'*X (boosts speed for n>>px) <optional>	    
%#									    
%# OUTPUT:    B    px x m  Y_on_X regression coefficients		    
%#	      C    m  x A  Y loadings					    
%#	      P    px x A  X loadings					    
%#	      T    n  x A  X scores                        <standardized>   
%#	      U    n  x A  Y scores					    
%#	      R    px x A  X weights					    
%#	      R2X  1  x A  cum. % X-variance accounted for		    
%#	      R2Y  1  x A  cum. % Y-variance accounted for		    
%#									    
%# WARNING:   In some cases, the number of factors extracted from X will    
%#	      be smaller than A because there is no covariance left of y    
%#	      with X. The column size of the output variables C, P, T, U,   
%#	      R, R2X and R2Y is then smaller than A, and the obtained B     
%#	      corresponds to a less than the A-factors model. 	    	    
%#									    
%# AUTHOR:    Sijmen de Jong, 15/6/1997					    
%# 	      Unilever Research Laboratorium, Vlaardingen, The Netherlands 
%#	      Copyright (c) 1997 for ChemoAC				    
%#	      Dienst FABI, Vrije Universiteit Brussel			    
%#	      Laarbeeklaan 103, 1090-Brussel Jette, Belgium		    
%#									    
%# VERSION: 1.1 (28/02/1998)						    
%#									    
%# TEST:      Vita Centner						    

function [B,C,P,T,U,R,R2X,R2Y]=plssim(X,Y,A,S,XtX);

[n,px] = size(X); [n,m] = size(Y);   				                % size of the input data matrices
if nargin<5, S = []; end, if isempty(S), S=(Y'*X)'; end		% if XtX not inputted, S=[]; always when S=[] then S=(Y'*X)'
if nargin<4, XtX=[]; end					                             % if S is not inputted, XtX=[];
if isempty(XtX) & n>3*px, XtX = X'*X; end			          % when XtX=[] and X is very "tall", the booster XtX is calculated
if nargin<3, A=10; end, A = min([A px n-1]);			       % if A is not inputted, then the defaul A is min[10 px n-1]
T = zeros(n ,A); U = T;						                            % initialization of variables
R = zeros(px,A); P = R; V = R;
C = zeros(m ,A); 
R2Y = zeros(1,A);
z = zeros(m,1); v = zeros(px,1); 
if n>px, S0 = S; end, StS = S'*S;				                % SIMPLS algorithm
nm1 = n-1;
tol = 0;
for a = 1:A
  StS = StS-z*z'; 
  [Q,LAMBDA] = eig(StS); 
  [lambda,j] = max(diag(LAMBDA)); 
  q = Q(:,j(1));
  r = S*q;
  t = X*r;
  if isempty(XtX), p = (t'*X)'; else p = XtX*r; end
  if n>px, d = sqrt(r'*p/nm1); else d = sqrt(t'*t/nm1); end
  if d<tol,      
	disp(' ')
    disp('WARNING: the required number of factors (A) is too high !')
	disp('Less PLS factors were extracted from the data in the PLSSIM program !') 
	disp(' ')
% 	break,
  else tol=max(tol,d/1e8);
  end
  v = p-V(:,1:max(1,a-1))*(p'*V(:,1:max(1,a-1)))'; v = v/sqrt(v'*v); 
  z = (v'*S)'; 
  S = S-v*z'; 
								% save results
  V(:,a) = v;
  R(:,a) = r/d; 						% X weights
  P(:,a) = p/(d*nm1); 						% X loadings
  T(:,a) = t/d;							% X scores
  U(:,a) = Y*q;							% Y scores
  C(:,a) = q*(lambda(1)/(nm1*d)); 				% Y loadings
  R2Y(1,a) =  lambda(1)/d;					% Y-variance accounted for
end
clear StS V LAMBDA Q p q r t v z;
% if d<tol,
%  A=a-1; a=A; T=T(:,1:A); U=U(:,1:A); R=R(:,1:A); P=P(:,1:A); C=C(:,1:A);
% end
while a>1
  U(:,a) = U(:,a)-T(:,1:a-1)*(U(:,a)'*T(:,1:a-1)/nm1)'; 
  a=a-1; 
end
B = R*C';							% B-coefficients of the regression Y on X
if isempty(XtX), sumX2=sum(X.^2); else sumX2 = sum(diag(XtX)); end
R2X = 100*nm1/sum(sumX2)*cumsum(sum(P.^2)); 
R2Y = 100/nm1/sum(sum(Y.^2))*cumsum(R2Y(1:A).^2);




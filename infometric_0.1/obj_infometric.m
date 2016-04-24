function [f, g] = obj_infometric(L, sX, sY, tX, lambda)

% Input: 
% L: D*d transformation matrix
% sX: instance matrix for the source domain data. size: N*D
% sY: label vector for the source domain data. size: N*1
% tX: instance matrix for the target domain data. size: M*D
% lambda: regularization parameter
% 
% Output:
% f: objective function value
% g: gradient on L
% 
% yuanshi@usc.edu
% 2013/9/5
%

% mutual info in the target domain
% do "discriminative clustering" on the target domain
[mi, grad] = compute_mutual_infoL(L, sX, sY, sX);

% mutual info in both domains
% make it hard to distinguish the two domains
all_data = [sX; tX];
all_label = [ones( size(sX,1), 1); zeros( size(tX,1), 1)];  % domain label

if lambda > 0
    [mi2, grad2] = compute_mutual_infoL(L, all_data, all_label, all_data);
    
    f = -(mi - lambda*mi2);
    g = -(grad - lambda*grad2);
else
    f = -mi;
    g = -grad;
end
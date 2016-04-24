function [ err, pangles ] = subspace_error(A, B)
%SUBSPACE_ERROR      Subspace error between two subspaces. 
%
%usage 
%  [err, pangles]  = subspace_error(A, B) 
%
%input
%  A    Basis for subspace A (column vectors). 
%  B    Basis for subspace B (column vectors).    
% 
%output 
%  err      Subspace error (mean of sin^2 of the principle angles)
%  pangles  Principle angles between subsace A and B. 
%
%author
%  paul.buenau@tu-berlin.de 

assert(size(A,1) == size(B, 1), 'A and B must have the same number of rows');

A = orth(A);
B = orth(B);

pangles = svd(A'*B);

err = 1 - mean(pangles.^2);
pangles = acos(pangles);

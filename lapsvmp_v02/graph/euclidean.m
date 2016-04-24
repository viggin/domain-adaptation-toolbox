function D = euclidean(A,B)
% {euclidean} computes the Euclidean distance.
%
%      D = euclidean(A,B)
%      
%      A: M-by-P matrix of M P-dimensional vectors 
%      B: N-by-P matrix of M P-dimensional vectors
% 
%      D: M-by-N distance matrix
%
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it
%         * based on the code of Vikas Sindhwani, vikas.sindhwani@gmail.com

if (size(A,2) ~= size(B,2))
    error('A and B must be of same dimensionality.');
end

if (size(A,2) == 1) % if dim = 1...
    A = [A, zeros(size(A,1),1)];
    B = [B, zeros(size(B,1),1)];
end

aa=sum(A.*A,2);
bb=sum(B.*B,2);
ab=A*B';

D = real(sqrt(repmat(aa,[1 size(bb,1)]) + repmat(bb',[size(aa,1) 1]) -2*ab));

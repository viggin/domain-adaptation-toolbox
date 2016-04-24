function A = adjacency(options,X)
% {adjacency} computes the graph adjacency matrix.
%
%      A = adjacency(options,X)
%
%      options: a structure with the following fields
%               options.NN: number of nearest neighbors to use
%               options.GraphDistanceFunction: 'euclidean' | 'cosine'
%               options.GraphWeights: 'distance' | 'binary' | 'heat'
%               options.GraphWeightParam: width for 'heat' kernel
%                                         (if set to 0, it uses the mean 
%                                          edge length distance among
%                                          neighbors) 
%      X: N-by-D data matrix (N examples, D dimensions)
%
%      A: sparse symmetric N-by-N adjacency matrix
%
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it
%         * based on the code of Vikas Sindhwani, vikas.sindhwani@gmail.com

n=size(X,1);
p=2:(options.NN+1);

if size(X,1)<500 % block size: 500
    step=n;
else
	step=500;
end

idy=zeros(n*options.NN,1);
DI=zeros(n*options.NN,1);
t=0;
s=1;

for i1=1:step:n    
    t=t+1;
    i2=i1+step-1;
    if (i2>n) 
        i2=n;
    end

    Xblock=X(i1:i2,:);  
    dt=feval(options.GraphDistanceFunction,Xblock,X);
    [Z,I]=sort(dt,2);
	 	    
    Z=Z(:,p)'; % it picks the neighbors from 2nd to NN+1th
    I=I(:,p)'; % it picks the indices of neighbors from 2nd to NN+1th
    [g1,g2]=size(I);
    idy(s:s+g1*g2-1)=I(:);
    DI(s:s+g1*g2-1)=Z(:);
    s=s+g1*g2;
end 

I=repmat((1:n),[options.NN 1]);
I=I(:);

if strcmp(options.GraphDistanceFunction,'cosine') % only positive values
    DI=DI.*(DI<1);
end

switch options.GraphWeights
    
    case 'distance' 
        A=sparse(I,idy,DI,n,n);
    
    case 'binary'
        A=sparse(I,idy,1,n,n);

    case 'heat'
        if options.GraphWeightParam==0 % default (t=mean edge length)
            t=mean(DI(DI~=0)); % actually this computation should be
                               % made after symmetrizing the adjacecy
                               % matrix, but since it is a heuristic, this
                               % approach makes the code faster.
        else
            t=options.GraphWeightParam;
        end

        A=sparse(I,idy,exp(-DI.^2/(2*t*t)),n,n);    
                
    otherwise
        error('Unknown weight type');
end

A=A+((A~=A').*A'); % symmetrize

function K = calckernel(options,X1,X2)
% {calckernel} computes the Gram matrix of a specified kernel function.
% 
%      K = calckernel(options,X1)
%      K = calckernel(options,X1,X2)
%
%      options: a structure with the following fields
%               options.Kernel: 'linear' | 'poly' | 'rbf' 
%               options.KernelParam: specifies parameters for the kernel 
%                                    functions, i.e. degree for 'poly'; 
%                                    sigma for 'rbf'; can be ignored for 
%                                    linear kernel 
%      X1: N-by-D data matrix of N D-dimensional examples
%      X2: (it is optional) M-by-D data matrix of M D-dimensional examples
% 
%      K: N-by-N (if X2 is not specified) or M-by-N (if X2 is specified)
%         Gram matrix
%
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it
%         * based on the code of Vikas Sindhwani, vikas.sindhwani@gmail.com

kernel_type=options.Kernel;
kernel_param=options.KernelParam;

n1=size(X1,1);
if nargin>2
    n2=size(X2,1);
end

switch kernel_type

    case 'linear'
        if nargin>2
            K=X2*X1';
        else
            K=X1*X1';
        end

    case 'poly'
        if nargin>2
            K=(X2*X1').^kernel_param;
        else
            K=(X1*X1').^kernel_param;
        end

    case 'rbf'
        if nargin>2
            K = exp(-(repmat(sum(X1.*X1,2)',n2,1) + ...
                repmat(sum(X2.*X2,2),1,n1) - 2*X2*X1') ...
                /(2*kernel_param^2));
        else
            P=sum(X1.*X1,2);
            K = exp(-(repmat(P',n1,1) + repmat(P,1,n1) ...
                - 2*X1*X1')/(2*kernel_param^2));
        end

    otherwise
        error('Unknown kernel function.');
end


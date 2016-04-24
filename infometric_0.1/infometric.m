function L = infometric(L0, sX, sY, tX, lambda)

% Unsupervised domain adaptation
%
% Input: 
% L0: D*d transformation matrix
% sX: instance matrix for the source domain data. size: N*D
% sY: label vector for the source domain data. size: N*1
% tX: instance matrix for the target domain data. size: M*D
% lambda: regularization parameter
% 
% Output:
% L: D*d transformation matrix 
% 
% Reference:
% Shi, Y., & Sha, F. (2012). Information-Theoretical Learning of
% Discriminative Clusters for Unsupervised Domain Adaptation. ICML 2012.
% 
%
% Notice: for practice use, you need to properly preprocess the data and 
% tune lambda and d. For details, please refer to the paper
%
% If you have any question, please contact
% yuanshi@usc.edu
% 2013/9/5
%

L = L0;
d = size(L,2);

MAXITER = 100;
stepsize = 1;

f = zeros(1,MAXITER);

for iter = 1:MAXITER
    [f(iter), g] = obj_infometric(L, sX, sY, tX, lambda);
    
    fprintf('%g  %f  %g\n', iter, f(iter), stepsize);
    
    % adjust the stepsize adaptively
    if iter > 1
        if f(iter) < f(iter-1)
            stepsize = stepsize*1.1;
        else
            stepsize = stepsize*0.5;
        end
        
        if stepsize > 20
            stepsize = 20;
        end
    end
    
    % gradient descent
    L = L - stepsize*g;
    
    % trace constraint
    ratio = d / trace(L'*L);
    L = sqrt(ratio)*L;
    
    % stopping condition
    if iter > 10 && max(abs( diff( f(iter-3:iter) ) ) )< 1e-6*abs(f(iter))
        fprintf('no more progress..%g %f\n', iter, f(iter));
        break;
    end
end
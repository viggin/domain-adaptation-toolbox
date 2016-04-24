function classifier = lapsvmp(options,data)
% {lapsvmp} trains a Laplacian SVM classifier in the primal.
%     
%      classifier = lapsvmp(options,data)
%
%      options: a structure with the following fields
%               options.gamma_A: regularization parameter (ambient norm)
%               options.gamma_I: regularization parameter (intrinsic norm)
%
%               [optional fields]
%               options.Cg: {0,1} i.e. train with Newton's method or PCG
%                           (default=0)
%               options.MaxIter: maximum number of iterations (default=200)
%               options.Hinge: {0,1} i.e. train a LapSVM (1) or LapRLSC (0)
%                              (default=1)
%               options.UseBias: {0,1} i.e. use or not a bias (default=0)
%               options.InitAlpha: if it's 0, the initial weights are null;
%                                  if it's <0, they are randomly taken;
%                                  otherwise it is the initial vector
%                                  (default=0)
%               options.InitBias: the initial bias (default=0)
%               options.NewtonLineSearch: {0,1} i.e. use or not exact line
%                                         search with Newton's method
%                                         (default=0).
%               options.NewtonCholesky: {0,1} i.e. use or not Cholesky
%                                       factorization and rank 1 updates of
%                                       the Heassian (default=0 for LapSVM.
%                                       If set to 1, the  Hessian must be
%                                       positive definite). 
%               options.CgStopType: {0,1,2,3,4,5,6,7,-1}
%                                   the data-based stopping criterion of cg
%                                   only (default=0), where:
%                                   0: do not stop
%                                   1: stability stop
%                                   2: validation stop
%                                   3: stability & validation stop
%                                   4: gradient norm (normalized)
%                                   5: prec. gradient norm (normalized)
%                                   6: mixed gradient norm (normalized)
%                                   7: relative objective function change
%                                  -1: debug (it saves many stats)
%               options.CgStopIter: number of cg iters after which check
%                                   the stopping condition.
%               options.CgStopParam: the parameter for the selected
%                                    CgStopType. If CgStopType is:
%                                    0: ignored
%                                    1: percentage [0,1] of tolerated
%                                       different decisions between two
%                                       consecutive checks
%                                    2: percentage [0,1] of error rate
%                                       decrement required  between two
%                                       consecutive checks
%                                    3: the two params above (i.e. it is an
%                                       array of two elements)
%                                    4: the minimum gradient norm
%                                    5: the minimum prec. gradient norm 
%                                    6: the minimum mixed gradient norm
%                                    7: relative objective function change
%                                       between two consecutive checks
%                                    -1: ignored
%                                    see the code for default values.
%               options.Verbose: {0,1} (default=1)
%
%      data: a structure with the following fields
%            data.X: a N-by-D matrix of N D-dimensional training examples
%            data.K: a N-by-N kernel Gram matrix of N training examples
%            data.Y: a N-by-1 label vector in {-1,0,+1}, where 0=unlabeled
%                    (it is in {0,+1} in the case of One-Class SVM/LapSVM.
%            data.L: a N-by-N matrix of the Laplacian
%
%            [other fields]
%            data.Kv: a V-by-N kernel Gram matrix of V validation examples,
%                     required if CgStopType is 2 or 3 (validation check)
%            data.Yv: a V-by-1 vector of {-1,+1} labels for the validation
%                     examples, required if CgStopType is 2 or 3
%                     (validation check)
%
%      classifier: structure of the trained classifier (see the
%                  'saveclassfier' function). 
%
% Author: Stefano Melacci (2012)
%         mela@dii.unisi.it
%         * the One-Class extension is joint work with Salvatore Frandina,
%           salvatore.frandina@gmail.com
%         * the original code structure was based on the primal SVM code of 
%           Olivier Chapelle, olivier.chapelle@tuebingen.mpg.de 

n=length(data.Y);
nn=nnz(data.Y==-1);

% initializing option structure with default values
if ~isfield(options,'Verbose'),           options.Verbose=1; end
if ~isfield(options,'Hinge'),             options.Hinge=1; end
if ~isfield(options,'Cg'),                options.Cg=0; end
if ~isfield(options,'MaxIter'),           options.MaxIter=200; end
if ~isfield(options,'UseBias'),           options.UseBias=0; end
if ~isfield(options,'InitAlpha'),         options.InitAlpha=false; end
if ~isfield(options,'InitBias'),          options.InitBias=0; end
if ~isfield(options,'NewtonLineSearch'),  options.NewtonLineSearch=0; end
if ~isfield(options,'NewtonCholesky'),    options.NewtonCholesky=0; end
if ~isfield(options,'CgStopType'),        options.CgStopType=0; end
switch options.CgStopType    
    case 0 % none
        if ~isfield(options,'CgStopParam') || ...
            isempty(options.CgStopParam), options.CgStopParam=0; end
        if ~isfield(options,'CgStopIter'), 
            options.CgStopIter=options.MaxIter+1; end         
    case 1 % stability
        if ~isfield(options,'CgStopParam') || ...
            isempty(options.CgStopParam), options.CgStopParam=0.015; end
        if ~isfield(options,'CgStopIter') || isempty(options.CgStopIter), 
            options.CgStopIter=round(sqrt(n)/2); end
    case 2 % validation
        v=length(data.Yv);      
        if ~isfield(options,'CgStopParam') || ...
            isempty(options.CgStopParam), options.CgStopParam=1/v; end
        if ~isfield(options,'CgStopIter') || isempty(options.CgStopIter), 
            options.CgStopIter=round(sqrt(n)/2); end           
    case 3 % stability & validation
        v=length(data.Yv);     
        if ~isfield(options,'CgStopParam') || ...
            isempty(options.CgStopParam), options.CgStopParam=...
                                                           [1/v,0.015]; 
        end
        if ~isfield(options,'CgStopIter'),
            options.CgStopIter=round(sqrt(n)/2); end
    case 4 % gradient norm
        if ~isfield(options,'CgStopParam') || ...
            isempty(options.CgStopParam), options.CgStopParam=1e-8; end
        if ~isfield(options,'CgStopIter') || isempty(options.CgStopIter),
            options.CgStopIter=1; end
    case 5 % preconditioned gradient norm
        if ~isfield(options,'CgStopParam') || ...
            isempty(options.CgStopParam), options.CgStopParam=1e-8; end
        if ~isfield(options,'CgStopIter') || isempty(options.CgStopIter),
            options.CgStopIter=1; end
    case 6 % mixed gradient norm
        if ~isfield(options,'CgStopParam') || ...
            isempty(options.CgStopParam), options.CgStopParam=1e-8; end
        if ~isfield(options,'CgStopIter') || isempty(options.CgStopIter),
            options.CgStopIter=1; end
    case 7 % relative objective function decrease
        if ~isfield(options,'CgStopParam') || ...
            isempty(options.CgStopParam), options.CgStopParam=1e-6; end
        if ~isfield(options,'CgStopIter') || isempty(options.CgStopIter),
            options.CgStopIter=1; end  
    case -1 % debug
        if ~isfield(options,'CgStopParam') || ...
            isempty(options.CgStopParam), options.CgStopParam=0; end
        if ~isfield(options,'CgStopIter') || isempty(options.CgStopIter), 
            options.CgStopIter=options.MaxIter+1; end        
    otherwise
        error('Invalid CgStopType.');         
end
if nn==0, oc=1; else oc=0; end

% initial alpha vector
if length(options.InitAlpha)>1
    alpha=options.InitAlpha;
else
    if options.InitAlpha==0
        alpha=[];
    else
        alpha=randn(n,1);
        alpha=alpha./norm(alpha);
    end
end

% checking common error conditions
if oc==1 && ~options.UseBias 
    error('One-Class SVM requires the UseBias option to be turned on.');
end    
if options.UseBias && options.LaplacianNormalize
    error(['The current implementation does not support a normalized ' ...
           'Laplacian when the UseBias options is turned on.']);
end
    
% initial bias
b=options.InitBias;

switch options.Cg
    case 0
        % newton
        [alpha,b,t,sec,lsiters]=newton(options,data,alpha,b,oc);
        stats=[];
    case 1
        % pcg
        [alpha,b,t,sec,lsiters,stats]=pcg(options,data,alpha,b,oc);
    otherwise
        error('Invalid solver specified in the field .Cg');
end

svs=find(alpha~=0);
classifier=saveclassifier('lapsvmp',svs,alpha(svs), ...
                          data.X(svs,:),b,options,sec,t,lsiters,stats);

function [alpha,b,t,sec,lsiters] = newton(options,data,alpha,b,oc)
% {newton} trains the classifier using the Newton's method.

tic
n=length(data.Y);
labeled=data.Y~=0;
l=nnz(labeled);
gamma_A=options.gamma_A;
gamma_I=options.gamma_I;

% initial seeding
if isempty(alpha)
    alpha=zeros(n,1); 
    Kalpha=zeros(n,1);
else
    Kalpha=data.K*alpha;
end

t=0;
lr=0;
sv=false(n,1);
if nargout>4, lsiters=zeros(options.MaxIter,1); end
                   
if gamma_I~=0, LK=data.L*data.K; end

while 1  
    
    if options.Hinge        
        sv_prev=sv;         
        hloss=sparse([],[],[],n,1,l);
        hloss(labeled)=1-data.Y(labeled).*(Kalpha(labeled,:)+b);    
        sv=hloss>0;
        nsv=nnz(sv);
    else
        sv_prev=sv; 
        sv=labeled;
        nsv=l; 
    end

    if options.Verbose
        if ~options.Hinge
            hloss=sparse([],[],[],n,1,l);
            hloss(labeled)=1-data.Y(labeled).*(Kalpha(labeled,:)+b);
        end
        if gamma_I~=0
            obj=(gamma_A*alpha'*Kalpha+sum(hloss(sv).^2)+...
                gamma_I*Kalpha'*data.L*Kalpha+oc*b)/2;
        else
            obj=(gamma_A*alpha'*Kalpha+sum(hloss(sv).^2)+oc*b)/2;
        end
 
        fprintf('[t=%d] obj=%f nev=%d lr=%.4f\n', [t full(obj) nsv lr]);
    end

    % goal conditions
    if t>=options.MaxIter, break, end           
    if isequal(sv_prev,sv), break, end 
    
    t=t+1;

    IsvK=sparse([],[],[],n,n,nsv*n);
    IsvK(sv,:)=data.K(sv,:);
        
    % computing new alphas
    onev=ones(1,n);

    if gamma_I==0 % SVM (sparse solution)
        if options.UseBias
            alpha_new=zeros(n,1);
            alpha_b_new=[0,onev(1:nsv);onev(1:nsv)',...
                         gamma_A*speye(nsv)+IsvK(sv,sv)]\ ...
                         [oc/(2*gamma_A);data.Y(sv)];
            alpha_new(sv)=alpha_b_new(2:end);
            b_new=alpha_b_new(1);
        else
            alpha_new=zeros(n,1);
            alpha_new(sv)=(gamma_A*speye(nsv)+IsvK(sv,sv))\data.Y(sv);
            b_new=0;
        end

    else % LapSVM
        
        % inversion by factorization
        if options.NewtonCholesky
            if t==1
                % compute the Cholesky factorization of the Hessian
                if options.UseBias
                    sumKsv=sum(data.K(sv,:));
                    hess=chol([nsv,sumKsv;sumKsv',...
                               data.K*(gamma_A*speye(n)+IsvK+gamma_I*LK)]);
                    alpha_b_new=hess\(hess'\([sum(data.Y(sv))-oc/2; ...
                                              data.K(:,sv)*data.Y(sv)]));
                    alpha_new=alpha_b_new(2:end);
                    b_new=alpha_b_new(1);
                else
                    hess=chol(data.K*(gamma_A*speye(n)+IsvK+gamma_I*LK));
                    alpha_new=hess\(hess'\(data.K(:,sv)*data.Y(sv)));
                    b_new=0;
                end
                LK=[];
            else
                % update the Cholesky factorization of the Hessian
                sv_diff=~(sv&sv_prev);
                sv_add=find(sv&sv_diff)';
                sv_rem=find(sv_prev&sv_diff)';
                if options.UseBias
                    if ~isempty(sv_add)
                        for i=sv_add
                            hess=cholupdate(hess,[1;data.K(:,i)],'+');
                        end
                    end
                    if ~isempty(sv_rem)
                        for i=sv_rem
                            hess=cholupdate(hess,[1;data.K(:,i)],'-');
                        end
                    end
                    alpha_b_new=hess\(hess'\([sum(data.Y(sv))-oc/2; ...
                                              data.K(:,sv)*data.Y(sv)]));
                    alpha_new=alpha_b_new(2:end);
                    b_new=alpha_b_new(1);
                else
                    if ~isempty(sv_add)
                        for i=sv_add
                            hess=cholupdate(hess,data.K(:,i),'+');
                        end
                    end
                    if ~isempty(sv_rem)
                        for i=sv_rem
                            hess=cholupdate(hess,data.K(:,i),'-');
                        end
                    end
                    alpha_new=hess\(hess'\(data.K(:,sv)*data.Y(sv)));
                    b_new=0;
                end
            end
        else         
            % inversion without factorization
            IsvY=sparse([],[],[],n,1,nsv); 
            IsvY(sv)=data.Y(sv);             
            if options.UseBias
                alpha_b_new=([0,onev;sv,...
                             gamma_A*speye(n)+IsvK+gamma_I*LK])\ ...
                             [oc/(2*gamma_A);IsvY];
                alpha_new=alpha_b_new(2:end);
                b_new=alpha_b_new(1);             
            else                
                alpha_new=(gamma_A*speye(n)+IsvK+gamma_I*LK)\IsvY;
                b_new=0;
            end 
            
        end
    end

    % step
    if options.NewtonLineSearch && (options.Hinge || nnz(alpha)>0)
        step=alpha_new-alpha;
        step_b=b_new-b;
        
        [lr,Kalpha,lsi]=linesearch(data,labeled,step,step_b,Kalpha,b,...
                                   gamma_A,gamma_I,[],[],options.Hinge,...
                                   oc);        
        alpha=alpha+lr*step;
        b=b+lr*step_b;
        lsiters(t)=lsi;
    else
        alpha=alpha_new;
        b=b_new; 
        lr=1;
        lsiters(t)=0;   
        Kalpha=data.K*alpha;     
    end
end

if nargout>4, lsiters=lsiters(1:t); end
sec=toc;


function [alpha,b,t,sec,lsiters,stats] = pcg(options,data,alpha,b,oc)
% {pcg} trains the classifier using preconditioned conjugate gradient.

tic
n=length(data.Y);
labeled=data.Y~=0;
unlabeled=data.Y==0;
l=nnz(labeled);
u=n-l;

if isfield(data,'Yv')
    v=length(data.Yv);
end
gamma_A=options.gamma_A;
gamma_I=options.gamma_I;

if isempty(alpha)
    alpha=zeros(n,1);
    Kalpha=zeros(n,1);
    if gamma_I~=0, LKalpha=zeros(n,1); else LKalpha=[]; end
    go=data.Y-b*labeled;
    obj0=(sum((1-data.Y(labeled)*b*options.UseBias).^2)+oc*b)/2;
    if options.UseBias, go_b=sum(data.Y(labeled)-b)-oc/2;
    else go_b=0; b=0; end
else
    Kalpha=data.K*alpha;
    if gamma_I~=0, LKalpha=data.L*Kalpha; else LKalpha=[]; end
    out=sparse([],[],[],n,1,l);
    out(labeled)=Kalpha(labeled)+b;
    if options.Hinge    
        sv=false(n,1);
        sv(labeled)=(data.Y(labeled).*out(labeled)<1);
    else
        sv=labeled;
    end
    go=-gamma_A*alpha;
    if gamma_I~=0, go=go-gamma_I*LKalpha; end
    obj0=(Kalpha'*(-g0)+sum((out(sv)-data.Y(sv)).^2)+oc*b)/2;
    go(sv)=go(sv)-(out(sv)-data.Y(sv));
    
    if options.UseBias, go_b=sum(data.Y(sv)-Kalpha(sv)-b)-oc/2;
    else go_b=0; b=0; end
    
end
d=go; % initial search direction
d_b=go_b;
Kgo=data.K*go;
Kstep=Kgo;

t=0;
stats=[];

switch options.CgStopType    
    case 4, ng0=sqrt(sum(Kgo.^2)+go_b^2);
    case 5, ngp0=sqrt(sum(go.^2)+go_b^2);
    case 6, ngm0=sqrt(sum(Kgo'*go)+go_b^2);        
    case -1
        stats=zeros(options.MaxIter+1,5+n+1);
        ng0=sqrt(sum(Kgo.^2)+go_b^2);
        ngp0=sqrt(sum(go.^2)+go_b^2);
        ngm0=sqrt(sum(Kgo'*go)+go_b^2);
        stats(t+1,:)=[t,obj0/obj0,ng0/ng0,ngp0/ngp0,ngm0/ngm0,alpha',b];
end

if nargout>3, lsiters=zeros(options.MaxIter,1); end
valerr_prev=1;
obj_prev=obj0;
yfx_unlabeled_prev=false(u,1);

while 1
    t=t+1;

    % goal condition: maximum number of iterations
    if t>options.MaxIter, t=t-1; break, end
    
    % do an exact line search   
    [lr,Kalpha,lsi,LKalpha]=linesearch(data,labeled,d,d_b,Kalpha,b,...
                                       gamma_A,gamma_I,Kstep,LKalpha,...
                                       options.Hinge,oc);   
                                   
    % goal condition: converged to optimal solution
    if lr==0, t=t-1; break, end
    
    alpha=alpha+lr*d;
    b=b+lr*d_b;
    lsiters(t)=lsi;

    % compute new precgradient and objective   
    out=sparse([],[],[],n,1,l);
    out(labeled)=Kalpha(labeled)+b;
    
    if options.Hinge    
        sv=false(n,1);
        sv(labeled)=(data.Y(labeled).*out(labeled)<1);
    else
        sv=labeled;
    end    
        
    g=gamma_A*alpha;
    if gamma_I~=0, g=g+gamma_I*LKalpha; end
    
    if options.Verbose || options.CgStopType==7 || options.CgStopType==-1
        obj=(Kalpha'*(g)+sum((out(sv)-data.Y(sv)).^2)+oc*b)/2;
        if options.Verbose        
            fprintf('[t=%d] obj=%f nev=%d lr=%.4f\n', ...
                     [t full(obj) nnz(sv) lr]);          
        end
    end
        
    g(sv)=g(sv)+(out(sv)-data.Y(sv));
    
    if options.UseBias, g_b=sum(out(sv)-data.Y(sv))+oc/2; else g_b=0; end
    
    gn=-g;
    gn_b=-g_b;

    % goal condition: data based
    if mod(t-1,options.CgStopIter)==(options.CgStopIter-1)
        switch options.CgStopType
            case 1 % stability
                yfx_unlabeled=sign(Kalpha(unlabeled)+b);
                diff_unlabeled=1-nnz(yfx_unlabeled==yfx_unlabeled_prev)/u;
                if diff_unlabeled<options.CgStopParam, break, end
                yfx_unlabeled_prev=yfx_unlabeled;
            case 2 % validation
                valerr=1-nnz(sign(data.Kv*alpha+b)==data.Yv)/v;   
                if (valerr>(valerr_prev-options.CgStopParam)), break, end
                valerr_prev=valerr;
            case 3 % stability & validation
                valerr=1-nnz(sign(data.Kv*alpha+b)==data.Yv)/v;    
                yfx_unlabeled=sign(Kalpha(unlabeled)+b);
                diff_unlabeled=1-nnz(yfx_unlabeled==yfx_unlabeled_prev)/u;
                if (valerr>(valerr_prev-options.CgStopParam(1))) && ...
                    diff_unlabeled<options.CgStopParam(2), break, end
                valerr_prev=valerr;
                yfx_unlabeled_prev=yfx_unlabeled;
            case 4 % gradient norm
                ng=sqrt(sum(gn.^2)+gn_b^2);
                if ng/ng0<options.CgStopParam, break, end           
            case 5 % preconditioned gradient norm
                ngp=sqrt(sum(go.^2)+go_b^2);
                if ngp/ngp0<options.CgStopParam, break, end                
            case 6 % mixed gradient norm
                ngm=sqrt(Kgo'*go+go_b^2);
                if ngm/ngm0<options.CgStopParam, break, end    
            case 7 % relative objective function decrease
                if (obj_prev-obj)<options.CgStopParam, break, end
                obj_prev=obj;
        end
    end
    
    Kgn=data.K*gn; % multiply by the preconditioner
    
    % debug
    if options.CgStopType==-1
        ng=sqrt(sum(Kgn.^2)+gn_b^2);
        ngp=sqrt(sum(gn.^2)+gn_b^2);
        ngm=sqrt(Kgn'*gn+gn_b^2);
        stats(t+1,:)=[t,obj/obj0,ng/ng0,ngp/ngp0,ngm/ngm0,alpha',b];
    end
 
    % Polack-Ribiere update with automatic restart
    be=max(0,(Kgn'*(gn-go)+gn_b*(gn_b-go_b))/(Kgo'*go+go_b^2));
    
    d=be*d+gn;
    d_b=be*d_b+gn_b;
    
    Kstep=be*Kstep+Kgn;

    go=gn;
    go_b=gn_b;
    Kgo=Kgn;
end
sec=toc;
if nargout>4, lsiters=lsiters(1:t);
    if nargout>5 && ~isempty(stats), stats=stats(1:t,:); end
end


function [lr,Kalpha,lsi,LKalpha] = linesearch(data,labeled,step,step_b, ...
                                              Kalpha,b,gamma_A,gamma_I, ...
                                              Kstep,LKalpha,hinge,oc)
% {linesearch} does a line search in direction step.

act=step~=0; % the set of points for which alpha change (active)
if isempty(Kstep)
    Kstep=data.K(:,act)*step(act);
end

% precomputations
stepKstep=step(act)'*Kstep(act); 
stepKalpha=step(act)'*Kalpha(act);
if gamma_I~=0
    KstepL=Kstep'*data.L;
    KstepLKalpha=KstepL*Kalpha;
    KstepLKstep=KstepL*Kstep;
end
     
out=Kalpha(labeled)+b;
outstep=Kstep(labeled)+step_b;
out_minus_Y=out-data.Y(labeled);

% breakpoints
if hinge
    sv=(1-out.*data.Y(labeled))>0;
    deltas=-out_minus_Y./outstep;
    deltas(deltas<0)=0;
    [deltas,deltas_map]=sort(deltas);
    lab=length(deltas);
    i=find(deltas>0,1);
    lsi=1;
else
    sv=true(length(out),1); 
    lsi=1;
end

% intercepts
if gamma_I~=0
    left=outstep(sv)'*out_minus_Y(sv)+gamma_A*stepKalpha+...
         gamma_I*KstepLKalpha + (oc/2)*step_b;
    right=outstep(sv)'*outstep(sv)+gamma_A*stepKstep+...
          gamma_I*KstepLKstep;
else
    left=outstep(sv)'*out_minus_Y(sv)+gamma_A*stepKalpha + (oc/2)*step_b;
    right=outstep(sv)'*outstep(sv)+gamma_A*stepKstep;    
end

% first minimum
zcross=-left/right;
if right<=0, zcross=0; end

if hinge && ~isempty(i)
    if sv(deltas_map(i))==true
        if 0<=zcross && zcross<deltas(i), not_got_it=0;
        else not_got_it=1; end
    else
        if 0<=zcross && zcross<=deltas(i), not_got_it=0;
        else not_got_it=1; end
    end

    while not_got_it
        % updating support vectors
        j=-2*sv(deltas_map(i))+1;
        sv(deltas_map(i))=~sv(deltas_map(i));

        % updating intercepts
        left=left+j*outstep(deltas_map(i))*out_minus_Y(deltas_map(i));
        right=right+j*outstep(deltas_map(i))*outstep(deltas_map(i)); 

        % computing minimum
        zcross=-left/right;

        % goal conditions
        if i==lab, break, end    
        if sv(deltas_map(i))==true
            if sv(deltas_map(i+1))==true
                if deltas(i)<zcross && zcross<deltas(i+1), 
                    not_got_it=0; end
            else
                if deltas(i)<zcross && zcross<=deltas(i+1), 
                    not_got_it=0; end
            end
        else
            if sv(deltas_map(i+1))==true
                if deltas(i)<=zcross && zcross<deltas(i+1), 
                    not_got_it=0; end
            else
                if deltas(i)<=zcross && zcross<=deltas(i+1), 
                    not_got_it=0; end
            end
        end

        i=i+1;
        lsi=lsi+1;
    end
end
lr=zcross;

if lr<0, lr=0; return, end % converged
Kalpha=Kalpha+lr*Kstep;
if nargout>3 && gamma_I~=0, LKalpha=LKalpha+lr*KstepL'; end

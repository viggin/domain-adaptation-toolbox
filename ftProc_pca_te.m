function [Xnew] = ftProc_pca_te(model,X)
%ftProc_pca_te Extract the PCA feature based on the trained coefficients.
%See FTPROC_PCA_TR.
%	[Xnew] = ftProc_pca_te(model,X)

%	Copyright 2015 Ke YAN, Tsinghua Univ. http://yanke23.tk, xjed09@gmail.com

Xnew = [];
if ~isempty(X)
	Xnew = bsxfun(@minus,X,model.mu);
	Xnew = Xnew * model.W_prj;
end

end
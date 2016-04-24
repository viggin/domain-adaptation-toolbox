function [Xnew] = ftProc_kpca_te(model,X)
%ftProc_kpca_te Extract the KPCA feature based on the trained coefficients.
%See FTPROC_KPCA_TR.
%	[Xnew] = ftProc_kpca_te(model,X)

%	Copyright 2015 Ke YAN, Tsinghua Univ. http://yanke23.tk, xjed09@gmail.com

Xnew = [];
if ~isempty(X)
	K = model.kerFun(X,model.Xtrain);
	Xnew = K * model.W_prj;
end

end
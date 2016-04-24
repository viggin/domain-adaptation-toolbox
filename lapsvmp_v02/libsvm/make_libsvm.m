% {make_libsvm} compiles Libsvm and the mex interface for Libsvm (C++).
%     
%      make_libsvm
%
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it

delete *.obj
delete *.o
mex -c svm.cpp
if exist('svm.obj','file')>0
    mex mexGramSVMTrain.cpp svm.obj
elseif exist('svm.o','file')>0
    mex mexGramSVMTrain.cpp svm.o
end

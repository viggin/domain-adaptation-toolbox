function Y = svmlfwd(net, X, Y)
% SVMLFWD - Wrapper for SVMlight: Prediction
%   
%   YPRED = SVMLFWD(NET, X)
%   Predict labels YPRED for points X using the SVMlight stored in NET.
%   Each row of X is one data point, X is of size [N D] for N points of
%   dimensionality D. The result Y is a column vector of size [N 1]. NET
%   is a structure holding information on the trained SVMlight model, as
%   it is output by SVMLTRAIN. In particular, the SVMlight model stored
%   in file NET.fname is used to classify the points X.
%   YPRED = SVMLFWD(NET, X, Y) also provides the actual target for the
%   test data. When using this syntax, SVMlight will print out the
%   accuracy on the test set, as well as precision and recall.
%   SVMLFWD(NET, FNAME), where FNAME is the name of an SVMlight data
%   file, computes the predictions for test data that are stored in file
%   FNAME.
%
%   See also SVML, SVMLTRAIN, SVM_CLASSIFY
%

% 
% Copyright (c) by Anton Schwaighofer (2002)
% $Revision: 1.2 $ $Date: 2002/02/19 12:23:16 $
% mailto:anton.schwaighofer@gmx.net
% 
% This program is released unter the GNU General Public License.
% 

error(nargchk(2, 3, nargin));
error(consist(net, 'svml'));

if nargin<3,
  Y = [];
end
fname = net.fname;
if ischar(X),
  testdata = X;
  deleteData = 0;
else
  testdata = [fname '.testdata'];
  svmlwrite(testdata, X, Y);
  deleteData = 1;
end
preddata = [fname '.preddata'];
status = svm_classify(net.options, testdata, net.fname, preddata);
Y = svmlread(preddata);
delete(preddata);
if deleteData,
  delete(testdata);
end
if status~=0,
  error(sprintf('Error when calling SVMlight. Status = %i', status));
end

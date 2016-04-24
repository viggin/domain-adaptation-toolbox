function demsvml1()
% DEMSVML1 - Demo program for Matlab SVMlight wrapper
%   
%   This program shows how Thorsten Joachims SVMlight test examples can
%   be executed using the Matlab SVMlight wrapper, plus some additional
%   features.
%   In order to execute this demo, SVMlight must be installed (on the
%   execution path or in the current directory) and the two
%   subdirectories 'example1' and 'example2' must be present in the
%   current directory.
%   
%   Please check the source code for more info.
%
%   See also SVML, SVMLTRAIN, SVMLFWD, SVMLWRITE, SVMLREAD, SVM_LEARN,
%            SVM_CLASSIFY
%

% 
% Copyright (c) by Anton Schwaighofer (2002)
% $Revision: 1.3 $ $Date: 2002/08/09 20:27:33 $
% mailto:anton.schwaighofer@gmx.net
% 
% This program is released unter the GNU General Public License.
% 

fprintf('In order to execute this demo, SVMlight must be installed (on the\n');
fprintf('execution path or in the current directory) and the two\n');
fprintf('subdirectories ''example1'' and ''example2'' must be present in the\n');
fprintf('current directory.\n');
fprintf('Have a look at the source code for details...\n');
fprintf('\nPress the Return key to continue, or Control-C to cancel.\n');

% We do basically the same as Joachims example2, just that we load
% the data into Matlab beforehand and use the SVML structure
[Y, X] = svmlread('example1/train.dat');
net = svml('example1/model1', 'Kernel', 0, 'C', 1);
net = svmltrain(net, X, Y);
% Also read out Joachims' test data
[Ytest, Xtest] = svmlread('example1/test.dat');
% Compute prediction on the test data
Xtest = Xtest(end-10:end,:);
Ypred = svmlfwd(net, Xtest);
% We can also get SVMlight to print out accuracy and precision/recall, by
% providing the test targets as well:
Ypred = svmlfwd(net, Xtest, Ytest);
% 
% 
% % The SVMlight wrapper can also directly access data files in SVMlight's
% % format, both during training and testing. This is of course more
% % efficient than reading the files into Matlab and writing it back to
% % disk again, if the data are already available in SVMlights format.
% net = svml('example1/model2', 'Kernel', 0, 'C', 1);
% % directly access training data
% net = svmltrain(net, 'example1/train.dat');
% % Compute prediction on the test data
% Ytest = svmlfwd(net, 'example1/test.dat');


% Now for the transductive SVM:
net = svml('example2/model1', 'Kernel', 0, 'C', 1, 'TransLabelFile', ...
           'example2/trans_predictions');
% Transductive learner is automatically invoked
% You may just as well load the data into Matlab variables
% [Y, X] = svmlread('example2/train_transduction.dat');
% and go net = svmltrain(net, X, Y);
net = svmltrain(net, 'example2/train_transduction.dat');
% Read out the labels for the unclassified training examples
transLabels = svmlread('example2/trans_predictions');


% We now use only the simple calling routines SVM_LEARN and
% SVM_CLASSIFY. For the options, we only assume that SVMlight is
% installed in the current directory '.'

% Joachims' first example (inductive SVM):
options = svmlopt('ExecPath', '.');
svm_learn(options, 'example1/train.dat', 'example1/model');
svm_classify(options, 'example1/test.dat', 'example1/model', ...
             'example1/predictions');
% The file example1/model holds the trained SVM model,
% example1/predictions are the labels for the test points

% Joachims' second example (transductive SVM)
svm_learn(options, 'example2/train_transduction.dat', 'example2/model');
svm_classify(options, 'example2/test.dat', 'example2/model', ...
             'example2/predictions');


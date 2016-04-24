function classifier = rlsc(options,data)
% {rlsc} trains a Regularized Least Square Classifier (RLSC).
%     
%      classifier = rlsc(options,data)
%
%      Note: This function is only a wrapper to 'lapsvmp', where the use of
%            a Hinge function for labeled data is disabled and the
%            contribute of unlabeled data is discarded. See 'lapsvmp'
%            for the description of all required parameters. All parameters
%            related to unlabeled data are obviously not checked.
%
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it

options.Hinge=0;
options.gamma_I=0;
classifier=lapsvmp(options,data);
classifier.name='rlsc';

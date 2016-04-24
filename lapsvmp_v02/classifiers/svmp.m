function classifier = svmp(options,data)
% {svmp} trains a SVM classifier on the primal.
%     
%      classifier = svmp(options,data)
%
%      Note: This function is only a wrapper to 'lapsvmp', where the
%            contribute of unlabeled data is discarded. See 'lapsvmp'
%            for the description of all required parameters. All parameters
%            related to unlabeled data are obviously not checked.
%
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it

options.Hinge=1;
options.gamma_I=0;
classifier=lapsvmp(options,data);
classifier.name='svmp';
function classifier = laprlsc(options,data)
% {laprlsc} trains a Laplacian RLSC.
%     
%      classifier = laprlsc(options,data)
%
%      Note: This function is only a wrapper to 'lapsvmp', where the use of
%            a Hinge function for labeled data is disabled. See 'lapsvmp'
%            for the description of all required parameters.
%
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it

options.Hinge=0;
classifier=lapsvmp(options,data);
classifier.name='laprlsc';
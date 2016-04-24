function classifier = oclapsvmp(options,data)
% {oclapsvmp} trains a One-Class Laplacian SVM classifier on the primal.
%     
%      classifier = oclapsvmp(options,data)
%
%      Note: This function is only a wrapper to 'lapsvmp', where the
%            contribute of unlabeled data is discarded. See 'lapsvmp'
%            for the description of all required parameters. All parameters
%            related to unlabeled data are obviously not checked.
%
% Author: Stefano Melacci (2012), Salvatore Frandina (2012)
%         mela@dii.unisi.it, salvatore.frandina@gmail.com

options.Hinge=1;
options.UseBias=1;
data.Y(data.Y==-1)=0;
classifier=lapsvmp(options,data);
classifier.name='oclapsvmp';
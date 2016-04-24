function [classifier] = saveclassifier(name,svs,alpha,xtrain,b,options,sec,t,lsiters,stats)
% {saveclassifier} generates structure with details of a kernel classifier.
%
%      classifier = saveclassifier(name,svs,alpha,xtrain,b,options)
%
%      classifier = saveclassifier(name,svs,alpha,xtrain,b,options,
%                                  sec,t,lsiters,stats)
%
%      name: string identifier of the classifier
%      svs: indices of the training points that are present on the
%           classifier function
%      alpha: coefficient vector (with |svs| elements)
%      xtrain: corresponding training examples
%      b: bias
%      options: the options data structure used to train the classifier,
%               that contains the regularization parameters, kernel 
%	            function, kernel parameters etc... 
%      
%      [optional]
%      sec: training time
%      t: number of iteration performed while training
%      lsi: vector with the number of line searches for each iteration
%      stats: matrix of misc statistics
%
%      classifier: is a data structure with the following fields
%                  classifier.name: name
%                  classifier.svs: svs
%                  classifier.alpha: alpha
%                  classifier.xtrain: xtrain
%                  classifier.b: b
%                  classifier.options: options
%
%                  [if optional arguments are present]
%                  classifier.iters: t
%                  classifier.lsiters: lsiters
%                  classifier.traintime: sec
%                  classifier.stats: stats
%
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it
%         * based on the code of Vikas Sindhwani, vikas.sindhwani@gmail.com 

classifier.name=name;
classifier.svs=svs;
classifier.alpha=alpha; 
classifier.b=b;
classifier.xtrain=xtrain;
classifier.options=options;
if exist('t','var') && ~isempty(t), classifier.iters=t; end
if exist('lsiters','var') && ~isempty(lsiters), classifier.lsiters=lsiters; end
if exist('sec','var') && ~isempty(sec), classifier.traintime=sec; end
if exist('stats','var') && ~isempty(stats), classifier.stats=stats; end

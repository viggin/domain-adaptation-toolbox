function resAll = doPredict(ft,target,cvObj)
% get the prediction accuracy. Masks for training and test data are in
% cvObj. The feature-level tranfer learning algorithm is in param.ftTrans.
% The logistic regression classifier is a part of the PRTools toolbox.

% masks
mas.tr = cvObj.training;  % data with labels
mas.te = cvObj.test; % data to predict
ftTr = ft(mas.tr,:);
ftTe = ft(mas.te,:);

pred1 = classf_('lr',ftTr,target(mas.tr),struct('lambda',1e0,'nIter',100),ftTe);
resAll = nnz(pred1 == target(mas.te))/length(pred1); % total acc

end

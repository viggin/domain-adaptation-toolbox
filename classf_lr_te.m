function [pred, prob] = classf_lr_te(model,Xtest)

Xtest = [ones(size(Xtest,1),1),Xtest];
prob = sigmoid(Xtest*model.thetas);
if size(model.thetas,2) == 1
	pred = prob >= .5;
	pred = 2-pred;
	prob = [prob,1-prob];
else
	[maxProb,pred] = max(prob,[],2);
end

end

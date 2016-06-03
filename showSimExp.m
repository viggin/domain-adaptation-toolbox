% synthetic dataset
function showSimExp

% close all
cvObj.NumTestSets = 1;
% rng('default');rng(0)

%% two discrete domains dataset
mu = {[0,0],[0,1];[1,-0.5]+randn(1,2)/10,[1,0.2]+randn(1,2)/10};
% mu = {[0,0],[0,1];[.5,.2],[.5,1.2]}; % less shift
Sigma = cell(size(mu)); 
for p = 1:numel(mu),
	X = diag(rand(2,1))*10;
	U = orth(rand(2));
	Sigma{p} = U'*X*U; 
end
% Sigma = {[1 .6; .6 1]*4};
% Sigma = repmat(Sigma,size(mu));
nSmps = 30;
nSmps = repmat(nSmps,size(mu));
[X,Y,domainFt,maLabeled] = genData(nSmps,mu,Sigma);

cvObj.training = maLabeled;cvObj.test = ~cvObj.training;
acc = doPredict(X,Y,cvObj);

r = 3; c = 3;
h = figure;
subplot(r,c,1)
draw1(X,Y,domainFt,{'x_1','x_2'},'Original',acc)
legend({'Pos. source data','Neg. source data','Pos. target data','Neg. target data'})

% TCA
param = []; param.kerName = 'lin';param.bSstca = 0;
param.mu = 1;param.m = 2;param.gamma = .1;param.lambda = 0;
[Xproj,transMdl] = ftTrans_tca(X,maLabeled,Y(maLabeled),maLabeled,param);
acc = doPredict(Xproj(:,1:2),Y,cvObj);
subplot(r,c,2)
draw1(Xproj,Y,domainFt,{'z_1','z_2'},'TCA',acc)

% MIDA
param = []; param.kerName = 'lin';param.kerSigma = 1e-1;param.bSup = 0;
param.mu = 1;param.m = 2;param.gamma = 1;
[Xproj,transMdl] = ftTrans_mida(X,domainFt,Y(maLabeled),maLabeled,param);
acc = doPredict(Xproj,Y,cvObj);
subplot(r,c,3)
draw1(Xproj,Y,domainFt,{'z_1','z_2'},'MIDA',acc)

% SA
param = []; param.pcaCoef = 2;
[Xproj,transMdl] = ftTrans_sa(X,maLabeled,Y(maLabeled),maLabeled,param);
acc = doPredict(Xproj,Y,cvObj);
subplot(r,c,4)
draw1(Xproj,Y,domainFt,{'z_1','z_2'},'SA',acc)

% ITL
param = []; param.pcaCoef = 1; param.lambda = 10;
[Xproj,transMdl] = ftTrans_itl(X,maLabeled,Y(maLabeled),maLabeled,param);
acc = doPredict(Xproj(:,1),Y,cvObj);
subplot(r,c,5)
draw1([Xproj,Xproj*0],Y,domainFt,{'z_1',''},'ITL',acc)

% GFK
param = []; param.dr = 1;
[Xproj,transMdl] = ftTrans_gfk(X,maLabeled,Y(maLabeled),maLabeled,param);
acc = doPredict(Xproj(:,1:2),Y,cvObj);
subplot(r,c,6)
draw1(Xproj,Y,domainFt,{'z_1','z_2'},'GFK',acc)

% PCA
param = []; param.pcaCoef = 2; param.kerName = 'lin';
[Xproj,transMdl] = ftTrans_pca(X,maLabeled,Y(maLabeled),maLabeled,param);
acc = doPredict(Xproj(:,1:2),Y,cvObj);
subplot(r,c,7)
draw1(Xproj,Y,domainFt,{'z_1','z_2'},'PCA',acc)

% Laplacian SVM
param = []; param.t = 0;
[pred,model,prob] = mdlTrans_lapsvm(X(maLabeled,:),Y(maLabeled),...
	X(~maLabeled,:),param);
acc = nnz(pred == Y(~maLabeled))/length(pred);

% Laplacian Ridge, use it for classification
param = []; param.t = 0;
[pred,model] = mdlTrans_lapridge(X(maLabeled,:),Y(maLabeled),...
	X(~maLabeled,:),param);
pred = (pred>1.5)+1;
acc1 = nnz(pred == Y(~maLabeled))/length(pred);
subplot(r,c,8)
title(sprintf('%s, %.2f%%\n%s, %.2f%%','Lap SVM',acc*100,'Lap Ridge',acc1*100))

set(gcf, 'position', [0 10 1350 900]);
% return
%}

%% continuous distribution change dataset
% if you need the SSA algorithm, please download at http://mloss.org/revision/download/851/
%{1
mu = {[0,0],[0 3]};
direc = {[0.3 1],[0.3 1];};
[nDomain,nCls] = size(mu);
Sigma = {[1 0; 0 2]/2};
Sigma = repmat(Sigma,size(mu));
nSmps = [100];
timeFt = (1:nSmps)'/10;
nSmps = repmat(nSmps,size(mu));
nSmpAll = sum(nSmps(:));

X = [];
Y = nan(nSmpAll,1);
domainFt = zeros(nSmpAll,nDomain*2);
maLabeled = [];
ord = Y;

for iDomain = 1:nDomain
	for iCls = 1:nCls
		mu1 = mu{iDomain,iCls};
		R1 = chol(Sigma{iDomain,iCls});
		z = repmat(mu1,nSmps(iDomain,iCls),1)+...
			timeFt*direc{iDomain,iCls} + randn(nSmps(iDomain,iCls),2)*R1;
		Y(size(X,1)+1:size(X,1)+nSmps(iDomain,iCls)) = iCls;
		domainFt(size(X,1)+1:size(X,1)+(nSmps(iDomain,iCls)),iDomain*2-1) = 1;
		domainFt(size(X,1)+1:size(X,1)+(nSmps(iDomain,iCls)),iDomain*2) = timeFt;
		ord((iDomain-1)*nDomain+iCls:nDomain*nCls:end) = size(X,1)+1:size(X,1)+(nSmps(iDomain,iCls));
		X = [X;z];
		maLabeled = [maLabeled;true(nSmps(iDomain,iCls)/2,1);false(nSmps(iDomain,iCls)/2,1)];
	end
end

X = X(ord,:);
Y = Y(ord,:);
domainFt = domainFt(ord,:);
maLabeled = maLabeled(ord,:);
maLabeled = maLabeled==1;

cvObj.training = maLabeled;cvObj.test = ~cvObj.training;
acc = doPredict(X,Y,cvObj);
figure,subplot(131)
draw2(X,Y,domainFt,{'x_1','x_2'},'Original',acc)
legend({'Pos. old data','Neg. old data','Pos. new data','Neg. new data'})

% SSA
param = []; param.d = 2;param.m = 1;
[Xproj,transMdl] = ftTrans_ssa(X,[],[],[],param);
acc = doPredict(Xproj(:,1),Y,cvObj);
subplot(132)
draw2(Xproj,Y,domainFt,{'z_1','z_2'},'SSA',acc)

% MIDA
param = []; param.kerName = 'lin';param.kerSigma = 1e-2;param.bSup = 0;
param.mu = 1;param.m = 2;param.gamma = 1;
param.doFtAug = 1;
[Xproj,transMdl] = ftTrans_mida(X,domainFt,Y(maLabeled),maLabeled,param);
acc = doPredict(Xproj,Y,cvObj);
subplot(133)
draw2(Xproj,Y,domainFt,{'z_1','z_2'},'MIDA',acc)
set(gcf, 'position', [0 300 1200 350]);
%}

%% back to the first datset, TSVM can be slow
% Transductive SVM
fprintf('running Transductive SVM...')
param = []; param.t = 0;
[pred,model,prob] = mdlTrans_tsvm(X(maLabeled,:),Y(maLabeled),...
	X(~maLabeled,:),param);
acc = nnz(pred == Y(~maLabeled))/length(pred);
figure(h);subplot(r,c,9)
title(sprintf('%s, %.2f%%','TSVM',acc*100))

end

function [X,Y,domainFt,maLabeled] = genData(nSmps,mu,Sigma)

nSmpAll = sum(nSmps(:));
maLabeled = false(sum(nSmps(:)),1);maLabeled(1:sum(nSmps(1,:))) = 1;
X = [];
Y = nan(nSmpAll,1);
[nDomain,nCls] = size(mu);
domainFt = false(nSmpAll,nDomain);
nDim = length(mu{1});

for iDomain = 1:nDomain
	for iCls = 1:nCls
		mu1 = mu{iDomain,iCls};
		R1 = chol(Sigma{iDomain,iCls});
		z = repmat(mu1,nSmps(iDomain,iCls),1) + randn(nSmps(iDomain,iCls),nDim)*R1/5;
		domainFt(size(X,1)+1:size(X,1)+nSmps(iDomain,iCls),iDomain) = 1;
		Y(size(X,1)+1:size(X,1)+nSmps(iDomain,iCls)) = iCls;
		X = [X;z];
	end
end
X = X-repmat(mean(X,1),size(X,1),1);

end

function draw1(X,Y,domainFt,axlb,name,acc)

% figure
hold on, %axis equal
co = {'r.','ro';'b+','bs';'g*','gd'};
nDomain = size(domainFt,2);
nCls = max(Y);
for iDomain = 1:nDomain
	for iCls = 1:nCls
		ma = domainFt(:,iDomain)&Y==iCls;
		plot(X(ma,1),X(ma,2),co{iDomain,iCls})
	end
end
xlabel(axlb{1})
ylabel(axlb{2})
title(sprintf('%s, %.2f%%',name,acc*100))

end

function draw2(X,Y,domainFt,axlb,name,acc)
% diff time smp

hold on
co = jet(nnz(domainFt(:,1)~=0&Y==1)); mo = '.o';
nDomain = size(domainFt,2)/2;
nCls = max(Y);

% for legend
for iDomain = 1:nDomain
	for iSmp = [1 90]
		for iCls = 1:nCls
			id = find(domainFt(:,iDomain*2-1)~=0&Y==iCls);
			[~,ord] = sort(domainFt(id,2),'ascend');
			ordr = 1:length(id); ordr(ord) = ordr;
			plot(X(id(iSmp),1),X(id(iSmp),2),mo(iCls),'color',co(ordr(iSmp),:))
		end
	end
end

for iDomain = 1:nDomain
	for iCls = 1:nCls
		id = find(domainFt(:,iDomain*2-1)~=0&Y==iCls);
		[~,ord] = sort(domainFt(id,2),'ascend');
		ordr = 1:length(id); ordr(ord) = ordr;
		for iSmp = 1:length(id)
			plot(X(id(iSmp),1),X(id(iSmp),2),mo(iCls),'color',co(ordr(iSmp),:))
		end
	end
end
xlabel(axlb{1})
ylabel(axlb{2})
title(sprintf('%s, %.2f%%',name,acc*100))

end

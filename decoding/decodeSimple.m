% this program comparing the decoding accuracy for random odor-Or
% interaction and optimized, maximum entropy encoding
% the input is just the Or response pattern, first use one layer direct
% decode the "compressed" signal

close all
clear

%% first set the basic parameters
nOdor = 50;
nRecp = 20;
sp = 3;
sig = 2;

NTrain = 1e5;
NTest = 2e4;

% learning parameters
beta = 1e-2;
MaxIter = 1e4;
MinBatch = 5e2;

trainData = trainDataGen(NTrain,nOdor,sp,sig);
testData = trainDataGen(NTest,nOdor,sp,sig);
% load the parameter for W with maximum entropy coding
dFolder = '/Users/shan/Dropbox/olfactionProject/data/gcmi/gcmi_Ndp';
file = 'gcmi_N50_R20_S3_sig2_2018-03-15.mat';
load(fullfile(dFolder,file));

% select a optimal w matrix
wInx = 8;               %select one w matrix
Wopt = reshape(allMat(:,wInx),[nRecp,nOdor]);

% random matrix
Wrnd = normrnd(mean(trainData(abs(trainData) > 0)),sig,nRecp,nOdor);

% get the input data
respOpt = exp(Wopt)*trainData./(1+exp(Wopt)*trainData);        % this is the input
respRnd = exp(Wrnd)*trainData./(1+exp(Wrnd)*trainData); 

%get the input of test
testOpt = exp(Wopt)*testData./(1+exp(Wopt)*testData); 
testRnd = exp(Wrnd)*testData./(1+exp(Wrnd)*testData); 

% using gradient descent to train the network weight of the decoder
% set the initial weight
w0 = sig*randn(nOdor*(nRecp+1),1);
% w0 = sig*randn(nOdor,nRecp);

%using Adam stocahstic gradient algorithm
% set minibatch
trainSize = size(trainData,2);
batchSize = 500;
numBatch = 1e3;
mnBatches = randi(trainSize,batchSize,numBatch);
cvnBatches = mat2cell(mnBatches, batchSize, ones(1, numBatch));
batchResp = @(nIter) respOpt(:,cvnBatches{mod(nIter,numBatch -1) + 1});
batchOdor = @(nIter) trainData(:,cvnBatches{mod(nIter,numBatch -1) + 1});

% define cost function
fhCost = @(w,nIter) targetFun(w,batchOdor(nIter),batchResp(nIter));

% options of Adam solver
sOpt = optimset('fmin_adam');
% sOpt.GradObj = 'off';
sOpt.TolX = 1e-4;
sOpt.MaxFunEvals = 1e4;
sOpt.Display = 'iter';
nEpochSize = 100;  
[xopt, fval, exitflag, output] = fmin_adam(fhCost,w0,0.01, [], [], [], nEpochSize, sOpt);

% test the performance on test set
wm = reshape(xopt,[nOdor,nRecp+1]);
Xhat = wm(:,1:end-1)*testOpt + wm(:,end)*ones(1,length(testOpt));
ftest = mean(sum((log(1+exp(Xhat)) - log(1+testData)).^2));


%% random measuring matrix situation

%train a random measuring matrix
batchResp = @(nIter) respRnd(:,cvnBatches{mod(nIter,numBatch -1) + 1});
batchOdor = @(nIter) trainData(:,cvnBatches{mod(nIter,numBatch -1) + 1});

% define cost function
fhCost = @(w,nIter) targetFun(w,batchOdor(nIter),batchResp(nIter));

% fhCost = @(w,nIter) targetFun(w,batchOdor(nIter),batchResp(nIter));
[xrnd, fval, exitflag, output] = fmin_adam(fhCost,w0,0.01, [], [], [], nEpochSize, sOpt);


wm = reshape(xopt,nOdor,nRecp+1);
Xhat = wm(:,1:end-1)*testOpt+ wm(:,end)*ones(1,length(testOpt));

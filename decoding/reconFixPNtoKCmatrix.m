% this program comparing the decoding accuracy for random odor-Or
% interaction and optimized, maximum entropy encoding
% the input is just the Or response pattern, first use one layer direct
% decode the "compressed" signal
% to incorporate the expansion from PN to KC, we limit the matrix from PN
% to KC fixed, but only alow the learning from KC to MBON

% this script is modified from "decodeSimple"

close all
clear

%% first set the basic parameters
nOdor = 50;
nRecp = 9;
sp = 2;
sig = 2;

NTrain = 1e5;
NTest = 2e4;

% learning parameters
beta = 1e-2;
MaxIter = 1e4;
MinBatch = 5e2;

% parameter for the PN to KC expansion
Nkc = 200;  % number of KC cells
rho = 0.1;  % sparsity of the expansition matrix

trainData = trainDataGen(NTrain,nOdor,sp,sig);
testData = trainDataGen(NTest,nOdor,sp,sig);
% load the parameter for W with maximum entropy coding
dFolder = '/Users/shan/Dropbox/olfactionProject/data/gcmi/gcmi_Ndp';
file = 'gcmi_N50_R9_S2_sig2_2018-03-20.mat';
load(fullfile(dFolder,file));

% select a optimal w matrix
wInx = 8;               %select one w matrix
Wopt = reshape(allMat(:,wInx),[nRecp,nOdor]);

% random matrix
Wrnd = normrnd(mean(trainData(abs(trainData) > 0)),sig,nRecp,nOdor);

% get the input data
respOpt = exp(Wopt)*trainData./(1+exp(Wopt)*trainData);        % this is the input
respRnd = exp(Wrnd)*trainData./(1+exp(Wrnd)*trainData); 

% the firt layer is a random sparse matrix
W1 = zeros(Nkc,nRecp);
W1(randperm(Nkc*nRecp,round(Nkc*nRecp*rho))) = rand(round(Nkc*nRecp*rho),1);
r_hat_opt = max(W1*respOpt-0.8,0);
r_hat_rnd = max(W1*respRnd-0.8,0);

%get the input of test
% testOpt = tanh(exp(Wopt)*testData./(1+exp(Wopt)*testData)); 
% testRnd = tanh(exp(Wrnd)*testData./(1+exp(Wrnd)*testData)); 

% using gradient descent to train the network weight of the decoder
% set the initial weight
% w0 = sig*randn(nOdor*(nRecp+1),1);
w0 = sig*randn(nOdor*(Nkc+1),1);

% w0 = sig*randn(nOdor,nRecp);

%using Adam stocahstic gradient algorithm
% set minibatch
trainSize = size(trainData,2);
batchSize = 1000;
numBatch = 1e3;
mnBatches = randi(trainSize,batchSize,numBatch);
cvnBatches = mat2cell(mnBatches, batchSize, ones(1, numBatch));
batchResp = @(nIter) r_hat_opt(:,cvnBatches{mod(nIter,numBatch -1) + 1});
batchOdor = @(nIter) trainData(:,cvnBatches{mod(nIter,numBatch -1) + 1});

% define cost function
fhCost = @(w,nIter) targetFun(w,batchOdor(nIter),batchResp(nIter));

% options of Adam solver
sOpt = optimset('fmin_adam');
% sOpt.GradObj = 'off';
sOpt.TolX = 1e-4;
sOpt.MaxFunEvals = 1e5;
sOpt.Display = 'iter';
sOpt.MaxIter = 1e4;
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

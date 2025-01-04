%% training a neural network to perform classification
% numTypes: the number of odor tyopes
% modified on 07/13/2017

function classifier(nOdor,nRecp,spar,sig,numTypes,varargin)

% parpool(2)
% parse input parameter, whether use random bumber seed
if nargin > 5
    rngSeed = varargin{1};
    rng(rngSeed)
else
    rng('shuffle')
end

% specify the input data type, random or cluster
inputType = 'cluster';


addpath(genpath('/lustre1/tangc_pkuhpc/ssqin/olfactoryCoding/optiInterMatrix/myFunc'))
spAll = 0.05:0.05:1;

% random number generator seeds

% ====== parameters of the network====
hiddenLayerSize = 100;  %set the hidden layer size
L = 1;                 %repeats
tsFun = 'tansig';      % transfer function, canbe pureliner, poslin,logsig,tansig
% ================================

% cycle through different number of hidden layers
allW1 = cell(length(spAll),L);
allW2 = cell(length(spAll),L);
allb1 = cell(length(spAll),L);
allb2 = cell(length(spAll),L);
allTestTrj = cell(length(spAll),L);
allTpr = cell(length(spAll),L);
allMisMatch = cell(length(spAll),L);
% ===========
% parameters on the measurement matrix, depends on input parameters
muW = -1;
sigW = 2;
% ============
for j0 = 1:L
    for i0 = 1:length(spAll)     
        sp = spAll(i0);  % input sparsity
        
        sparsity = spar;
        param.nRecep = nRecp;
        param.nOdor = nOdor;
        param.lSig = sig;
        param.lMu = 0;
        param.nSamp = 5e4;
        param.sparsity = true;
        param.spRatio = sparsity/param.nOdor;
        param.dType = 'lognorm';
        param.eig = [];
        param.corr = false;
        param.regularize = false;
        param.spType = 'absolute';          %type of sparisty, could be absolute or average
        param.noiseSig = 1e-2;              % introduce small noise on the response
        
        param.nPattern = 1e3;        % number of pattern, centroids
        param.withinPattSamp = 1e2;  % number of samples in each pattern
        param.patternStd = 0.1; 
        

        w = normrnd(muW,sigW,[nRecp,nOdor]);
        w_zero = rand(nRecp,nOdor);
        w(w_zero>sp) = -inf;
        
        if strcmp(inputType,'random')
%             trainData0 = zeros(param.nOdor,param.nSamp);
            eigVal = specifyEig(param.nOdor,param.eig);
            corrCoefMat = randCorrCoef('buildin',eigVal);
            trainData = genTrainData(param,corrCoefMat);
            NS = param.nSamp;
            targets = zeros(2,param.nSamp);  
            for k0 = 1:param.nSamp
                inx = find(trainData(:,k0) > 0);
                if all(inx <= param.nOdor/2)
                    targets(1,k0) = 1;
                elseif all(inx > param.nOdor/2)
                    targets(2,k0) = 1;
                else
                    if trainData(inx(1),k0) > trainData(inx(2),k0)
                    targets(1,k0) = 1;
                    else
                    targets(2,k0) = 1;
                    end
                end
              
            end
        elseif strcmp(inputType,'cluster')
            [trainData, targets] = genClusterInput(param,spar,numTypes);
            NS = param.nPattern*param.withinPattSamp;
        else
            error('the input data type can be only  random or cluster!')
        end
        
        % assign targets label, 2xNSamp matrix
%           
        
        
        % add some noise on the input o not
        
        if ~isempty(param.noiseSig)
            resp = reshape(exp(w),nRecp,nOdor)*trainData./(1+reshape(exp(w),nRecp,nOdor)*trainData)...
                + param.noiseSig*randn(param.nRecep,NS);
        else
            resp = reshape(exp(w),nRecp,nOdor)*trainData./(1+reshape(exp(w),nRecp,nOdor)*trainData);
        end
        
        % Create a Pattern Recognition Network
        
        net = patternnet(hiddenLayerSize);
        net.layers{1}.transferFcn = tsFun;
        
        net.trainParam.epochs = 2e3;
        net.trainParam.max_fail=20;    %increase validation check when fail


        % Set up Division of Data for Training, Validation, Testing
        net.divideParam.trainRatio = 70/100;
        net.divideParam.valRatio = 15/100;
        net.divideParam.testRatio = 15/100;

        % Train the Network
        [net,tr] = train(net,resp,targets,'useParallel','yes');
        
        allW1{i0,j0} = net.IW{1};
        allW2{i0,j0} = net.LW{2};
        allb1{i0,j0} = net.b{1};
        allb2{i0,j0} = net.b{2};
        allTestTrj{i0,j0} = tr.tperf;
        
        % Test the Network
%         outputs = net(resp);
%         errors = gsubtract(targets,outputs);
%         performance = perform(net,targets,outputs);

        % performanace on the test data
        tInd = tr.testInd;
        tstOutputs = net(resp(:,tInd));
        tstPerform = perform(net,targets(:,tInd),tstOutputs);
        
        allTpr{i0,j0} = tstPerform;
        
        % record all the ratio of mismatch
        [c,~,~,~] = confusion(targets(:,tInd),tstOutputs);
        allMisMatch{i0,j0} = c;

    end
end
% use a struct to save all the data
summData = struct('hiddenSize',hiddenLayerSize,'noiseSig',param.noiseSig,...
    'allTestTrj',allTestTrj,'allTpr',allTpr,'allW1',allW1,'allW2',allW2,...
    'allb1',allb1,'allb2',allb2,'tsFun',tsFun,'misMatch',allMisMatch);

% save the final trainging results
sName = ['classify_N',num2str(nOdor),'M',num2str(nRecp),'_HS',num2str(hiddenLayerSize),...
    '_sp',num2str(spar),'_numType',num2str(numTypes),'_',date,'.mat'];
save(sName,'summData');
end

function [trainData, targets] = genClusterInput(param,spar,numTypes)
% param      a struct specify all the parameters needed

        
NP = param.nPattern;       % number of pattern, centroids
NS = param.withinPattSamp;  % number of samples in each pattern
ds = param.patternStd;      % variation level within pattern


% set the centroid;
% trainCentroid = zeros(param.nOdor,NP);
eigVal = specifyEig(param.nOdor,param.eig);
corrCoefMat = randCorrCoef('buildin',eigVal);
trainCentroid = genTrainData(param,corrCoefMat);

% randomly assign label to these centroid
centroidLabel = zeros(numTypes,NP);
for i0 = 1:NP
    centroidLabel(randi(numTypes),i0) = 1;
end
% centroidLabel = binornd(1,1/numTypes,[1,NP]);
% centroidLabel = [centroidLabel;~centroidLabel];

% sample and assign label within each pattern
trainData = zeros(param.nOdor,NP*NS);
targets = zeros(numTypes,NP*NS);
for i0 = 1:NP
%     inx = find(trainCentroid(:,i0) > 0);
    trainData(trainCentroid(:,i0) > 0,(i0-1)*NS+1:1:i0*NS) = (1+ds*randn(spar,NS)).*trainCentroid(trainCentroid(:,i0) > 0,i0);
    targets(:,(i0-1)*NS+1:1:i0*NS) = repmat(centroidLabel(:,i0),1,NS);
end
end
% this program try to use CMA-ES algorithm to search for the optimal
% odor-or interaction matrix
% last revised 06/26/2017

clc
clear
global trainData numOdor numRecep regul
%add path when run on linux cluster
% addpath(genpath('/lustre1/tangc_pkuhpc/ssqin/olfactoryCoding/optiInterMatrix'))

%% this part specify training data type, optimization method and other stuff
sparsity = true;     %sparse data
dataType = 'lognorm'; %wide or narrow distribution
regul = false;         %if use regularization of entropy
correlation = false;   %if introduce correlation  in training data
%optMethod = 'wholeMat'; %optimize whole matrix, or just the distribution
optMethod = 'distr'; %optimize whole matrix, or just the distribution

%% set parameters
numSamp = 1e3;   
if sparsity
    numOdor = 50;  %number of odors
    numRecep = 5;  %number of receptors
    spars= 5; %average number of ligand in mixture
    trainData = zeros(numOdor,numSamp);   
else
    % non-sparse data
    numOdor = 20;  %number of odors
    numRecep = 4;  %number of receptors
    trainData = zeros(numOdor,numSamp);
    spars = numRecep;
end

% data type
if strcmp(dataType,'gauss')
       mu = 5;
       sig = 1;
%        trainData(inx) = abs(normrnd(mu,sig,[numSamp*sparsity,1]));    
       OdorMeanCon = mu;
       if correlation
           eigVal = specifyEig(numOdor,numRecep,[7,5,4,2]);
           corrCoefMat = randCorrCoef('buildin',eigVal);
           trainData = generateTrainData(numSamp,numOdor,corrCoefMat,mu,sig,'gauss',spars);
       else
           trainData = abs(normrnd(mu,sig,[numSamp,numOdor]));
       end
elseif strcmp(dataType,'lognorm')
        mu = 0;
        sig = 2;
        OdorMeanCon = exp(mu + sig^2/2); 
        % with correlation, although the method should be examined
        if correlation
            %eigVal = specifyEig(numOdor,numRecep,[7,5,4,3]);  % fixed the first 4 dimension for repeatable
            eigVal = specifyEig(numOdor,numRecep);  % fixed the first 4 dimension for repeatable
            corrCoefMat = randCorrCoef('buildin',eigVal);
            trainData = generateTrainData(numSamp,numOdor,corrCoefMat,mu,sig,spars);
        % default, without correlation
        else
            inx = datasample(1:numOdor*numSamp,numSamp*spars,'Replace',false);
            trainData(inx) = exp(normrnd(mu,sig,[numSamp*spars,1]));   
        end
end


%% optimization methods
%initial value
mu_w = -log(OdorMeanCon*spars); 
iniSig = 10; %inital standard deviation

if strcmp(optMethod,'wholeMat')  %optimize the whole matrix
    opts.StopFitness = 1e-2;
    opts.MaxIter = 1e4; %maximum number of iteration
    opts.LBounds = -15;
    opts.UBounds = 15;
%     opts.PopSize  = 100;
    
    w0 = normrnd(mu_w,sig,[numRecep*numOdor,1]);
    [wmin,fmin] = cmaes('entropyCost',w0,iniSig,opts);    
else %optimize distribution
    iniSig = [10;10;1]; %inital standard deviation
    opts.StopFitness = 1e-2;
    opts.MaxIter = 1e4; %maximum number of iteration
    opts.LBounds = [-15;0;0];
	pts.UBounds = [15;20;1];
    opts.PopSize  = 100;
    
    w0 = [mu_w,sig,0.1];
    [wmin,fmin] = cmaes('entropyCostDist',w0,iniSig,opts);    
end

% save optimal matrix and target function
%saveFile = ['out_',num2str(id),'.mat'];
save('out.mat','wmin','fmin')
%clear defined persistent variabels
clear entropyCost


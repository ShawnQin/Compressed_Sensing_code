function [wmin,fmin] = optWparam(numOdor,numRecp,sparsity,sig,rho,varargin)
% Search for the optimal mu_w and sig_w when fixing the sparsity of W
% modified from 'optMatrixCMA_v2_edited.m'
% N            dimension of input, or number of odor
% M            dimension of measurement, or number of receptors
% S            sparsity of input, an integer
% sig          standard deviation of input data
% varargin{1}  rng random number seed

% this program integrates different algorithms together, so that we can
% choose different method systematically

%% orgnize the input arguments
p = inputParser;
p.addRequired('numOdor',@isnumeric);
p.addRequired('numRecp',@isnumeric);
p.addRequired('sparsity',@isnumeric);
p.addRequired('sig',@isnumeric);
p.addRequired('rho',@isnumeric);

% p.addRequired('overlap',@isnumeric);
p.addParameter('numSamp',2e4,@isnumeric);  % maximum iteration times
p.addParameter('maxIter',1e1,@isnumeric);  % maximum iteration times
p.addParameter('popSize',100,@isnumeric);  % population size
p.addParameter('noiseSig',1e-3,@isnumeric);%noise level
p.addParameter('overlap',1e-3,@isnumeric); % degree of over lap
p.addParameter('lambda',2,@isnumeric);     % strength of panelty
p.addParameter('r0',0.2,@isnumeric);       % strength of panelty
p.addParameter('frac',0.5,@isnumeric);     % the fraction of inhibitory receptors
p.addOptional('rngSeed',@isnumeric);       % the rng seed

parse(p,numOdor,numRecp,sparsity,sig,rho,varargin{:})

numOdor = p.Results.numOdor;
numRecp = p.Results.numRecp;
numSamp = p.Results.numSamp;
sparsity = p.Results.sparsity;
rhow = p.Results.rho;

PopSize = max(4*(4+ceil(3*log(numOdor*numRecp))),100);
maxIter = p.Results.maxIter;
noiseSig = p.Results.noiseSig;
% overLap = p.Results.overlap;
rngSeed = p.Results.rngSeed;   %specify the random number seed from external fle


global param trainData
% global param trainData

%add information tool box when runinig on CLS cluster
% addpath(genpath('/lustre1/tangc_pkuhpc/lqy/olfaction'))
addpath(genpath('/lustre1/tangc_pkuhpc/ssqin/olfactoryCoding/optiInterMatrix/myFunc'))
%% ----------- global parameter -------------------------

param.nOdor = numOdor;            %number of total odors
param.nRecep = numRecp;           %number of receptors

param.dType = 'exponential';          %odor concentratio distribution, power law, lognormal, exponential
param.method = 'distr';             %optimization method, 'wholeMat', which use KDP estimator of entropy
                                    % 'gcmi' use gmci and KDP to estimate differential entropy
                                    % 'twobymul', use direct integration algorithm
param.nSamp = numSamp;              %number of odor mixture

param.lMu = 1;                      % mean for lognorm
param.lSig= sig;                    % std for lognorm
param.rho = rhow;                   % sparsity of W

param.sparsity = true;              %sparse input
param.corr = false;                 %correlation across odor concentration
param.regul = false;                %if use further reuglarization in cost function
param.noise = false;                %if add noise
param.noiseSig = noiseSig;          %specify the sigma of gaussian noise
param.eig = [numOdor/2,numOdor/4,numOdor/8,numOdor/16];                     %specify the first 4 component in data,column vectors
param.spRatio = sparsity/numOdor;   %sparse ratio, percentage of non-zero elements
param.h = 1;                        %Hill coefficient of response function
param.regularize = true;            %if apply regularization on W
param.constraint = true;            % introudce structural constraints on odorants

param.spType = 'absolute';          %type of sparisty, could be absolute or average

%% random number generator seed, revised on 04/05/2018
if isnumeric(rngSeed)
    rng(rngSeed)
else
    rng('shuffle')
end

%% options of cma-es alogrithm
opts = cmaes;
OdorMeanCon = exp(param.lMu(1) + param.lSig(1)^2/2);
mu_w = -log(OdorMeanCon*param.spRatio*param.nOdor);

eigVal = specifyEig(param.nOdor,param.eig);
corrCoefMat = randCorrCoef('buildin',eigVal);
trainData = genTrainData(param,corrCoefMat);


%% ------------ different optimization methods -------------------------%%
% default use modified KDP

if strcmp(param.method,'distr')  %only search for specific distribution
    
%     iniSig = 10;
    opts.MaxIter = maxIter;  %maximum number of iteration
    opts.LBounds = [-5*param.lSig;0.1*param.lSig];
	opts.UBounds = [5*param.lSig;5*param.lSig];
    
    w0 = [mu_w,param.lSig];
    [wmin,fmin,~,~,~,bestever] = cmaes('entrCostParamW',w0,[],opts);
   
else
    error('Method has to be either wholeMat or distr or anal or inhibit !')
end

%% ------------- save the results ---------------------- 
dFolder = './'; %current directory, when run on cluster

fStr = [num2str(param.nOdor),'x',num2str(param.nRecep),...
      '_indep_',param.method,'_sp',num2str(round(param.spRatio*param.nOdor)),'.mat'];
fName = fullfile(dFolder,fStr);
save(fName,'wmin','fmin','bestever')

end
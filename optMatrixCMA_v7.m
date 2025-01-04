function [wmin,fmin] = optMatrixCMA_v7(numOdor,numRecp,sparsity,sig,Cons,r0,varargin)
% This is the 7nd version of optMatrixCMA.m, modified from
% "optMatrixCMA_v2_edited.m"
% numOdor           dimension of input, or number of odor, 100
% numRecp           dimension of measurement, or number of receptors, 30
% sparsity          sparsity of input, number of odorants, an integer, 2
% sig               standard deviation of input data, lognorm dist., 2
% Cons              the contraint of input, used when adding correlation
% r0                relative baselne activity, default 0
% varargin{1}       number of sampling, default is 50


%% orgnize the input arguments
p = inputParser;
p.addRequired('numOdor',@isnumeric);
p.addRequired('numRecp',@isnumeric);
p.addRequired('sparsity',@isnumeric);
p.addRequired('sig',@isnumeric);
p.addRequired('Cons',@ischar); 
p.addRequired('r0',@isnumeric);
p.addParameter('numSamp',2e4,@isnumeric);   % total number of odor samples
p.addParameter('maxIter',2e3,@isnumeric);   % maximum iteration times
p.addParameter('popSize',100,@isnumeric);   % population size
p.addParameter('noiseSig',1e-3,@isnumeric); %noise level
p.addOptional('rngSeed',@isnumeric);        %the rng seed

parse(p,numOdor,numRecp,sparsity,sig,Cons,r0,varargin{:})

numOdor = p.Results.numOdor;
numRecp = p.Results.numRecp;
numSamp = p.Results.numSamp;
sparsity = p.Results.sparsity;
r0 = p.Results.r0;
PopSize = p.Results.popSize;
maxIter = p.Results.maxIter;
noiseSig = p.Results.noiseSig;
constraint = p.Results.Cons;
rngSeed = p.Results.rngSeed;

global param X trainData

%add information tool box when runinig on CLS cluster
% addpath(genpath('/lustre1/tangc_pkuhpc/ssqin/olfactoryCoding/optiInterMatrix/myFunc'))

%% ----------- global parameter -------------------------

param.nOdor = numOdor;            %number of total odors
param.nRecep = numRecp;           %number of receptors

param.dType = 'lognorm';     %odor concentratio distribution
param.method = 'wholeMat';   %optimization method, whole matrix or distribution, "wholeMat/distr"
param.function = false;
param.nSamp = numSamp;       %number of odor mixture

param.lMu = 0;               % mean for lognorm, default 0
param.lSig= sig;             % std for lognorm
% param.gMu = 5;               %mean for gaussian data
% param.gSig = 1;              %std for gaussian data

param.sparsity = true;       %sparse input
param.corr = false;          %correlation across odor concentration
%param.regul = true;         %if use further reuglarization in cost function
param.noise = false;         %if adding noise, for numeric stability
param.noiseSig = noiseSig;   %specify the sigma of gaussian noise
param.eig = [];              %specify the first 4 component in data,column vectors
param.spRatio = sparsity/param.nOdor;       %sparse ratio, percentage of non-zero elements
param.r_0 = r0;                 %alpha value = Rmax/R0 - 1 for inhibitory response
param.target = 'a';
param.regularize = true;      % regularize the weight, see the target function, entropy calculation for detail
param.spType = 'absolute';    % type of sparisty, could be absolute or average


% Determine the different contraints of odor ditribution
if strcmp(constraint,'none')
    param.constraint = false;

elseif strcmp(constraint,'group')
    param.constraint = true;
    cons = 0.2*rand(param.nOdor,param.nOdor);
    groupNum = 3;
    group = [floor(param.nOdor/groupNum).*ones(1,groupNum-1),param.nOdor-(groupNum-1)*floor(param.nOdor/groupNum)];
    matr = cell(1,groupNum);
    for i0 = 1:groupNum
        matr{i0} = 0.5+0.3*rand(group(i0),group(i0));
    end
    corre = blkdiag(matr{:});
    param.cos = (cons + corre)'-tril(cons + corre)';
    param.cos = param.cos(param.cos~=0)';
else
    param.constraint = true;
    if strcmp(constraint,'inde')
        alpha = 1;
        beta = 10;
    elseif strcmp(constraint,'corr')
        alpha = 10;
        beta = 1;
    elseif strcmp(constraint,'med')
        alpha = 2;
        beta = 2;
    elseif strcmp(constraint,'doub')
        alpha = 0.5;
        beta = 0.5;    
    else
        error('Constraint has to be "corr" or "med" or "doub" or "inde"!');
    end
    cons = betarnd(alpha,beta,[1,param.nOdor*(param.nOdor-1)/2]);
    param.cos = cons;
end

%% options of cma-es alogrithm
opts = cmaes;
OdorMeanCon = exp(param.lMu(1) + param.lSig(1)^2/2);
mu_w = -log(OdorMeanCon*param.spRatio*param.nOdor);
X = abs(normrnd(0,1e-10,param.nRecep,param.nSamp));

eigVal = specifyEig(param.nOdor,param.eig);
corrCoefMat = randCorrCoef('buildin',eigVal);
trainData = genTrainData(param,corrCoefMat);
if isnumeric(rngSeed)
	rng(rngSeed)
else
	rng('default')
end
%trainData = (ones(1e4,1)*[1,1,0])';
%trainData(trainData==0) = abs(normrnd(0,1e-10,1,length(find(trainData==0))));
%%------------------------------------------%%
 % a very small random variable to avoid response to be exactly 1 or 0

%% ------------ different optimization methods -------------------------%%
if strcmp(param.method,'wholeMat')
    % options
    %addpath('./gmciDistrN100M20S2sig2_new');
%     x = load(strcat('gcmi_distr_N100_R20_S2_sig',num2str(param.lSig),'_2018-08-28.mat'),'allParam');
%     x = x.allParam(:,1);
    iniSig = 20;     % std of inital w ditribution
    opts.StopFitness = [];
    opts.MaxIter = maxIter;   %maximum number of iteration
    opts.LBounds = [];
    opts.UBounds = [];
    opts.PopSize = PopSize;
%     nonZeroInx = randperm(param.nRecep*param.nOdor,round(param.nRecep*param.nOdor*x(3)));
%     w(nonZeroInx) = normrnd(x(1),abs(x(2)),[length(nonZeroInx),1]);
    %w0 = normrnd(mu_w,param.lSig,[param.nRecep*param.nOdor,1]);
    w0 = 6*rand(param.nRecep*param.nOdor,1)*param.lSig-3*param.lSig;
    [wmin,fmin,~,~,~,bestever] = cmaes('gcmiCost',w0,iniSig,opts);    
elseif strcmp(param.method,'distr')
%     iniSig = [10;10;1]; %inital standard deviation
    iniSig = [];
%     opts.StopFitness = 1e-3;
    opts.MaxIter = 1e4; %maximum number of iteration
    opts.LBounds = [-20;0;0];
	opts.UBounds = [5;3*param.lSig;1];
%    opts.PopSize  = 100;
    spar = 0.5;
    w0 = [0,param.lSig,spar];
%     w0 = [mu_w,param.lSig];
    [wmin,fmin,~,~,~,bestever] = cmaes('entropyCostDist',w0,iniSig,opts);
elseif strcmp(param.method,'twobymul')
    % only two odorants
    iniSig = 10;
    opts.StopFitness = [];
    opts.MaxIter = 2e3; %maximum number of iteration
    opts.LBounds = -100;
    opts.UBounds = [];
    opts.PopSize = 100;
     w0 = normrnd(mu_w,param.lSig,[param.nRecep*param.nOdor,1]);
    [wmin,fmin,~,~,~,bestever] = cmaes('entropycostAnav2',w0,iniSig,opts); 
elseif strcmp(param.method,'onetomul')
    % only one odorant
    iniSig = 10;
    opts.StopFitness = 0;
    opts.MaxIter = 1e3; %maximum number of iteration
    opts.LBounds = [];
    opts.UBounds = [];
    opts.PopSize = 50;
     w0 = normrnd(mu_w,param.lSig,[param.nRecep*param.nOdor,1]);
    [wmin,fmin,~,~,~,bestever] = cmaes('entropycostAnav1',w0,iniSig,opts); 
elseif strcmp(param.method,'inhibit') %consider the inhibitory response, not working too well
    opts.StopFitness = [];
    opts.MaxIter = 5e3; %maximum number of iteration for searching the Kij
    opts.LBounds =1e-4;
    opts.UBounds =1-1e-4;
    opts.PopSize = 100;
    iniSig =0.4;
    w0 = rand([param.nRecep*param.nOdor]);
    [wmin,fmin,~,~,~,bestever] = cmaes('entropyCostInhibit',w0,iniSig,opts); 
    wmin = reshape(wmin,param.nRecep*param.nOdor,1);
    r_0 = param.r_0;
    alpha = (1-r_0)/r_0;
    Sign = 2*(wmin>r_0)-1;
    wmin = ((1./wmin-1)/alpha).^(-Sign)./1e-3-1;
    wmin = [reshape(wmin,param.nRecep*param.nOdor,1),reshape(Sign,param.nRecep*param.nOdor,1)];
else
    error('Method has to be either wholeMat or distr or anal or inhibit !')
end

%% ------------- save the results ---------------------- 
dFolder = './'; %current directory, when run on cluster
if param.corr
    fStr = [num2str(param.nOdor),'x',num2str(param.nRecep),...
        '_corr_',param.method,'_sp',num2str(round(param.spRatio*param.nOdor)),'.mat'];
else
    fStr = [num2str(param.nOdor),'x',num2str(param.nRecep),...
        '_indep_',param.method,'_sp',num2str(round(param.spRatio*param.nOdor)),'.mat'];
end
fName = fullfile(dFolder,fStr);
save(fName,'wmin','fmin','bestever')
end

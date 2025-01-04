% this program comparing the classification accuracy for random odor-Or
% interaction and optimized, maximum entropy encoding
% the input is just the Or response pattern and related tag (odor valence), 
% to incorporate the expansion from PN to KC, we limit the matrix from PN
% to KC fixed, but only alow the learning from KC to MBON
% we also use global inhibition at the KC level to generate sparse
% representation

% Use simple Hebbian classifier
% this script is modified from "classifyRndPNtoKC"
% last revised on 08/26/2018

function HebbianClassifier(nOdor,nRecp,spar,sig,varargin)
% nOdor         number of odorants, an integer
% nRecp         number of receptors, an integer
% sig           odorant concentration, defaut 2
% numTypes      groups of lables, an interger larger than 1
% varargin{1}   extral input, such as the random number seed


% parse input parameter, whether use random number seed
if nargin > 5
    rngSeed = varargin{1};
    rng(rngSeed)
else
    rng('shuffle')
end

% specify the input data type, random or cluster
inputType = 'cluster';   % generating clustered odor clouds
classMethod = 'Hebbian';     % classification methods, can be logistic,svm,LDA,or multiClass
KCmethod = 'binary';       % activation function at KC level, can be relu, binary, linear

% add the search path, this is only useful when runing on the cluster
addpath(genpath('/lustre1/tangc_pkuhpc/ssqin/olfactoryCoding/optiInterMatrix/myFunc'))
spAll = 0.1:0.1:1;     %the range of sparsity of W    

% ====== parameters of the network======================
hiddenLayerSize = 1000;      %set the hidden layer size
L = 20;                     %repeats
% ======================================================


% =======================================================
% parameters on the measurement matrix, should be adjusted according to the
% number of odorant and ORN
% [muW,sigW,rhoW,fmin] = selectMeanSig(nOdor,spar);
muW = -1.55;
sigW = 1.85;
% ========================================================

% the parameters of the input
sparsity = spar;
param.nRecep = nRecp;
param.nOdor = nOdor;
param.lSig = sig;
param.lMu = 0;
% param.nSamp = 2e4;
param.sparsity = true;
param.spRatio = sparsity/param.nOdor;
param.dType = 'lognorm';
param.eig = [];
param.corr = false;
param.regularize = false;
param.spType = 'absolute';          %type of sparisty, could be absolute or average
param.noiseSig = 1e-6;              % introduce small noise on the response
        
param.nPattern = 200;        % number of pattern, centroids
param.withinPattSamp = 3;   % number of samples in each pattern
param.patternStd = 0.01;
param.sigw = 0.1;           %the distribution of PN-KC weight

spPoject = 6;             % average number of PN  each KC receive
KCsp = 0.1;               % KC response sparsity

% initialization
errRate = zeros(length(spAll),L);  % store the results
for j0 = 1:L
    for i0 = 1:length(spAll)     
        sp = spAll(i0);  % input sparsity
        
        w = normrnd(muW,sigW,[nRecp,nOdor]);
        w_zero = rand(nRecp,nOdor);
        w(w_zero>sp) = -inf;
        
        if strcmp(inputType,'cluster')
            [testSet, targets,centroid,centroidLabel] = genOdorClouds(param,spar);
            NS = param.nPattern*param.withinPattSamp;
            param.nSamp = NS;
        else
            error('the input data type can be only  random or cluster!')
        end
        
        % add some noise on the input or not
        % random truncated Gaussian noise
        pd = makedist('Normal');
        t = truncate(pd,0,inf);
        if ~isempty(param.noiseSig)
            % centrod response
            resp0 = reshape(exp(w),nRecp,nOdor)*centroid./(1+reshape(exp(w),nRecp,nOdor)*centroid)...
                + param.noiseSig*random(t,param.nRecep,param.nPattern);
            % testing response
            resp = reshape(exp(w),nRecp,nOdor)*testSet./(1+reshape(exp(w),nRecp,nOdor)*testSet)...
                + param.noiseSig*random(t,param.nRecep,NS);
        else
            resp0 = reshape(exp(w),nRecp,nOdor)*centroid./(1+reshape(exp(w),nRecp,nOdor)*centroid);
            resp = reshape(exp(w),nRecp,nOdor)*testSet./(1+reshape(exp(w),nRecp,nOdor)*testSet);
        end
        
        % define the matrix from PN to KC, a random projection matrix
        % response at the KC level,introduce global inhibition
        
%         W1 = zeros(hiddenLayerSize,nRecp);
        W1 = PNtoKC(nRecp,hiddenLayerSize,spPoject,'binary');
%         W1(randperm(hiddenLayerSize*nRecp,round(hiddenLayerSize*nRecp*spPoject))) ...
%             = normrnd(1,param.sigw,round(hiddenLayerSize*nRecp*spPoject),1);
%         W1(randperm(hiddenLayerSize*nRecp,round(hiddenLayerSize*nRecp*spPoject))) ...
%             = rand(round(hiddenLayerSize*nRecp*spPoject),1);

        % get the read out neuron weight using Hebbian supervised learning
        rhat0 = mean(resp0,2)/norm(mean(resp0,2));
        normResp0 = W1*(resp0 - rhat0*(rhat0'*resp0));  % from Abbott 2010 PNAS
        KCout0 = genKCresp(normResp0,KCmethod,1-KCsp);
        Wclass = (KCout0 - KCsp)*centroidLabel';
        
        
        % test the Hebbian classifier
        % first get the input
        rhat = mean(resp,2)/norm(mean(resp,2));
        normResp = W1*(resp - rhat*(rhat'*resp));  % from Abbott 2010 PNAS
%         normResp = (W1 - mean(W1(W1 > 0))*spPoject)*resp;  % Abbott 2017 Neuron
%         KCinput = W1*resp;
        
        % generate sparse representation at KC level
        KCout = genKCresp(normResp,KCmethod,1-KCsp);
        
        errRate(i0,j0) = mean(abs(targets' - sign((KCout' - KCsp)*Wclass)))/2;
        
    end
end

% plot the figure
% errorbar(spAll',mean(errRate,2),std(errRate,0,2))
% use a struct to save all the data
summData = struct('hiddenSize',hiddenLayerSize,'noiseSig',param.noiseSig,...
   'errorRate',errRate,'allSp',spAll,'KCsp',KCsp);

% save the final trainging results
sName = ['classify_',classMethod,'_N',num2str(nOdor),'M',num2str(nRecp),'_HS',...
    num2str(hiddenLayerSize),'_sp',num2str(spar),'_group',num2str(numTypes),'_',date,'.mat'];
save(sName,'summData');
end
%% generate KC response 
function KCout = genKCresp(resp,method,varargin)
% return the sparse KC output
% resp should be the total input to KC, linear function of PN/OSN
% method is a string, specified in "wta",'relu','linear'


% winner takes all
if strcmp(method,'binary')
    KCout = zeros(size(resp,1),size(resp,2));
    if nargin > 2      
        thd = varargin{1}; %threhold for the rectified linear output, 0 ~1
    else
        thd = 0.95;   % default only select the 5% strongest response
    end    
    
    
    tempThd = thd;
    for i0 = 1:size(resp,1)
        [F,X] = ecdf(resp(i0,:));
        inx = find(F>tempThd,1,'first');
        thd = X(inx);
        KCout(i0,resp(i0,:)>thd) = 1; 
    end
    
elseif strcmp(method,'relu')
    if nargin > 2      
        thd = varargin{1}; %threhold for the rectified linear output, 0 ~1
    else
        thd = 0.95;   % default only select the 5% strongest response
    end    
        % tune the threshold such that only 10% KC response
    [F,X] = ecdf(resp(:));
    inx = find(F>thd,1,'first');
    thd = X(inx);
    
%     MAXRESP = max(resp(:)) - thd;
%     KCout = max(0,resp - thd)/MAXRESP; %normalized to 0 ~ 1
    KCout = max(0,resp - thd);
        
elseif strcmp(method,'linear')
    KCout = resp;  %just return the original data
else
    error('method not support yet, has to be wta, relu, or linear!')
end
end

%% plot the sparsity dependent classification error
function plotFigure(summData)

% folder to save the figure
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% default graphic setting
defaultGraphicsSetttings

% color
Bu = brewermap(11,'Blues');    % blues

errorbar(summData.allSp',mean(summData.errorRate,2),std(summData.errorRate,0,2),...
    'o-','MarkerSize',12,'MarkerFaceColor',Bu(10,:),'Color',Bu(10,:),'LineWidth',2,...
    'CapSize',0)
lg = legend('N=100,M=20,sp=3,\sigma_c=2');
legend boxoff
set(lg,'FontSize',16)
xlabel('sparsity of W')
ylabel('classification error')

figNamePref = ['clasErr_PN_KCrand_N100M20_sp3_sig2_H500_diffSpW_logistic',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])
end

%% connection from PN to KC
function W = PNtoKC(nRecp,hiddenLayerSize,sp,weightType,varargin)
% nRecp                 number of receptors, PN
% hiddenLayerSize       number of KCs, an integer
% sp                    average number of PNs projected to KC
% weightType            could be gaussian, or uniform random, or binary
% varargin{1}           exact connection numbers each KC receive

% initialize the W
W = zeros(hiddenLayerSize,nRecp);
allConn = binornd(10,(sp-1)/11,[hiddenLayerSize,1]) + 1;
allConn(allConn>nRecp) = nRecp;


if strcmp(weightType,'binary')
    for i0 = 1:hiddenLayerSize
        W(i0,randperm(nRecp,allConn(i0))) = 1;
    end
elseif strcmp(weightType,'gaussian')
    % truncated Gaussian noise
    pd = makedist('Normal');
    t = truncate(pd,-2,2);
    for i0 = 1:hiddenLayerSize
        W(i0,randperm(nRecp,allConn(i0))) = (random(t,allConn(i0),1) + 2)/4;
    end
elseif strcmp(weightType,'uniform')
    for i0 = 1:hiddenLayerSize
        W(i0,:) = rand(allConn(i0),1);
    end
end
end

%% generate centroid odor mixture that are more uniform in the space
function [testSet, targets,centroid,centroidLabel] = genOdorClouds(param,spar,varargin)
% param      a struct specify all the parameters needed
% spar       sparsity of odor
% varargin{1}  spcecify the centroid concentration type, identical, or
% random

if nargin > 4
    cType = varargin{1};
else
    cType = 'random';   % default
end

        
NP = param.nPattern;        % number of pattern, centroids
NS = param.withinPattSamp;  % number of samples in each pattern
ds = param.patternStd;      % variation level within pattern


% set the centroid, which is chose uniformly cover the space

% step 1, get the index
inx = zeros(NP,spar);
for i0 = 1:NP
    inx(i0,:) = randperm(param.nOdor,spar);
end

% step 2, set the concentration, range -2sig ~ 2 sig
centroid = zeros(param.nOdor,NP);
if strcmp(cType,'identifcal')
    for i0 = 1:NP
        centroid(inx(i0,:),i0) = 1;
    end
elseif strcmp(cType,'random')
    for i0 = 1:NP
        centroid(inx(i0,:),i0) = exp(randn(1,spar)*param.lSig);
    end
    
end

% step 3 randomly assign label to these centroid, with -1 and 1
centroidLabel = rand(1,NP);
centroidLabel(centroidLabel>=0.5) = 1;
centroidLabel(centroidLabel<0.5) = -1;

% step 4, generate clouds around centroids
testSet = zeros(param.nOdor,NP*NS);
targets = zeros(1,NP*NS);
for i0 = 1:NP
    testSet(centroid(:,i0) > 0,(i0-1)*NS+1:1:i0*NS) = exp((1+ds*randn(spar,NS))).*centroid(centroid(:,i0) > 0,i0);
    targets((i0-1)*NS+1:1:i0*NS) = repmat(centroidLabel(i0),1,NS);
end
end 
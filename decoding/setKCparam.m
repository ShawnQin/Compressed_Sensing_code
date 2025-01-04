% this program try to set the optimal parameter of KC layer in our
% classification and reconstruction tasks
function setKCparam()

%% set the basica parameters
nOdor = 100;
nRecp= 20;
spar = 3;
sig = 2;
numTypes = 2;  % group the odors

spPoject = 7/nRecp;             %number of PN  each KC receive
% trainRatio = 0.8;
KCsp = 0.1;   
hiddenLayerSize = 500;  %set the hidden layer size
KCmethod = 'wta';



spW = 0.1:0.1:1;

%% the parameters for generating sparse input
% the parameters of the input
sparsity = spar;
param.nRecep = nRecp;
param.nOdor = nOdor;
param.lSig = sig;
param.lMu = 0;
param.nSamp = 1e4;
param.sparsity = true;
param.spRatio = sparsity/param.nOdor;
param.dType = 'lognorm';
param.eig = [];
param.corr = false;
param.regularize = false;
param.spType = 'absolute';          %type of sparisty, could be absolute or average
param.noiseSig = 1e-2;              % introduce small noise on the response
        
param.nPattern = 1e3;        % number of pattern, centroids
param.withinPattSamp = 10;  % number of samples in each pattern
param.patternStd = 0.1;
param.sigw = 0.1;      %the distribution of PN-KC weight

%%
% [muW,sigW,rhoW,fmin] = selectMeanSig(nOdor,spar);
muW = -2;
sigW = 2;
%% generate input data
[trainData, targets] = genClusterInput(param,spar,numTypes);
NS = param.nPattern*param.withinPattSamp;

%% response of ORN/PN
% resp = reshape(exp(w),nRecp,nOdor)*trainData./(1+reshape(exp(w),nRecp,nOdor)*trainData);
%% loop through W with different sparsity
ratioMissOdor = zeros(length(spW),1);
stdKC = zeros(length(spW),1);
ratioMissKC = zeros(length(spW),1);
allThd = zeros(length(spW),1);  % all the threshold
figure
hold on

for i0 = 1:length(spW)
    sp = spW(i0);  % input sparsity
        
    w = normrnd(muW,sigW,[nRecp,nOdor]);
    w_zero = rand(nRecp,nOdor);
    w(w_zero>sp) = -inf;
    resp = reshape(exp(w),nRecp,nOdor)*trainData./(1+reshape(exp(w),nRecp,nOdor)*trainData)...
        + param.noiseSig*randn(nRecp,param.nPattern*param.withinPattSamp);

    W1 = zeros(hiddenLayerSize,nRecp);
%     W1(randperm(hiddenLayerSize*nRecp,round(hiddenLayerSize*nRecp*spPoject))) ...
%             = normrnd(1,param.sigw,round(hiddenLayerSize*nRecp*spPoject),1);
    W1(randperm(hiddenLayerSize*nRecp,round(hiddenLayerSize*nRecp*spPoject))) ...
            = rand(round(hiddenLayerSize*nRecp*spPoject),1);
%     rhat = mean(resp,2)/norm(mean(resp,2));
    normResp = (W1 - mean(W1(W1 > 0))*spPoject)*resp;
%     normResp = W1*(resp - rhat*(rhat'*resp));  % from Abbott 2010 PNAS
%     totInput = W1*resp;
        
% generate sparse representation at KC level
    KCout = genKCresp(normResp,KCmethod,1-KCsp);
    z1 = sum(KCout > 0,1);
    z2 = sum(KCout > 0,2);
    stdKC(i0) = std(z1);
    ratioMissOdor(i0) = sum(z1==0)/size(trainData,2);
    ratioMissKC(i0) = sum(z2==0)/hiddenLayerSize;
    
    subplot(3,4,i0)
    histogram(z1)
end
hold off

figure
plot(spW',ratioMissOdor)
xlabel('sp of W')
ylabel('fraction missing odor')

figure
plot(spW',stdKC)
xlabel('sp of W')
ylabel('std of active KC')

figure
plot(spW',ratioMissKC)
xlabel('sp of W')
ylabel('fraction of useless KC')

end

%%
function KCout = genKCresp(resp,method,varargin)
% return the sparse KC output
% resp should be the total input to KC, linear function of PN/OSN
% method is a string, specified in "wta",'relu','linear'


% winner takes all
if strcmp(method,'wta')
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
    
    
%     num = varargin{1}; %precentage of each column
%     [~,Ix] = maxk(KCout,round(num*size(resp,1)),1);
%     for k0 = 1:size(KCout,2)
%         KCout(Ix(:,k0),k0) = 1; 
%     end
    
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
    
    MAXRESP = max(resp(:)) - thd;
    KCout = max(0,resp - thd)/MAXRESP; %normalized to 0 ~ 1
        
elseif strcmp(method,'linear')
    KCout = resp;  %just return the original data
else
    error('method not support yet, has to be wta, relu, or linear!')
end
end
%% generate input data function
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

%% select the optimal parameters for sparse matrix
function [muW,sigW,rhoW,fmin] = selectMeanSig(N,sp)
% this function return the optimal parameter of parameterized W

% load the data
dFoler = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';  % run local
% dFoler = '../';  % run on cluster
dName= fullfile(dFoler,'gcmi_distri_summData_16-Jul-2018.mat');
load(dName)
% allN = allN;
% allSp = allSp;
inx1 = find(allN ==N);
inx2 = find(allSp ==sp);

muW = meanW(inx1,inx2);
sigW = meanSigW(inx1,inx2);
rhoW = meanRho(inx1,inx2);
fmin = meanFmin(inx1,inx2);
end

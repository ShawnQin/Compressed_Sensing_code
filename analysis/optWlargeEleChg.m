% this program test how changing the dominant large matirx elements change
% the performance of optimal matrix
% In some case in our previous simulation, we observed there exist a largr
% bump in the sensitive part of W

close all
clear

%% first, specifify the parameters used in the simulation

param.nOdor = 50;            %number of total odors
param.nRecep = 20;           %number of receptors
sparsity = 3;                % number of odorant in a odor mixture
param.noiseSig = 1e-3;          %specify the sigma of gaussian noise
param.lSig= 2;                    % std for lognorm
param.nSamp = 2e4;              %number of odor mixture


param.dType = 'lognorm';          %odor concentratio distribution, power law, lognormal, exponential
param.method = 'gcmi';             %optimization method, 'wholeMat', which use KDP estimator of entropy
                                    % 'gcmi' use gmci and KDP to estimate differential entropy
                                    % 'twobymul', use direct integration algori
param.lMu = 0;                      % mean for lognorm
param.lambda = [];                 %the distribution used to generate correlation distance
param.c0 = 1;                       %minimum concentration, power law

param.sparsity = true;              %sparse input
param.corr = false;                 %correlation across odor concentration
param.noise = false;                %if add noise
param.eig = [];                     %specify the first 4 component in data,column vectors
param.spRatio = sparsity/param.nOdor;   %sparse ratio, percentage of non-zero elements
param.h = 1;                        %Hill coefficient of response function
param.regularize = true;            %if apply regularization on W
param.spType = 'absolute';          %type of sparisty, could be absolute or average


opts = cmaes;

%% generate training data and check its distribution
eigVal = specifyEig(param.nOdor,param.eig);
corrCoefMat = randCorrCoef('buildin',eigVal);
trainData = genTrainData(param,corrCoefMat);

% check the distribuiton, a histogram
figure
z = trainData(trainData> 0);
histogram(log(z))
xlabel('$\ln(w)$','interpreter','latex')
ylabel('count')

%% Second, load the optimal sensitivity matrix from data
inx = 1;   % default index of which matrix to test
[fileName, pathname] = uigetfile({'*.mat'},'File Selector');
fullName = fullfile(pathname,filesep,fileName);
load(fullName);

w0 = allMat(:,inx);
f0 = allfmin(inx);

% a histogram of the matrix
thr = -6;  % threshold of non-zero elements
ub = 3;  % upper bound, used in heatmap
figure
histogram(w0(w0>thr))

% a heat map of w
figure
w = reshape(w0,[param.nRecep,param.nOdor]);
imagesc(w,[thr,ub])

strongInx = w0 > ub;
% new matrix
wn = w0;
wn(strongInx) = mean(w0(w0>thr))+4;  % set the strong elements into average value

%% Third, calculate the entropy for optimal and new W
% using information tool box
mult = 1;
col = HShannon_KDP_initialization(mult);

% optimal matrix
resp = exp(w)*trainData./(1+exp(w)*trainData) + param.noiseSig*randn(param.nRecep,size(trainData,2));
MI = nonparanormal_info(resp');
H0 = 0;
for i0 = 1:param.nRecep
    H0 = H0 + HShannon_KDP_estimation(resp(i0,:),col);
end
H = -H0 + MI; 
disp(['the entropy for original optimal W is: ',num2str(H)])

% for the new matrix
w2 = reshape(wn,[param.nRecep,param.nOdor]);
resp = exp(w2)*trainData./(1+exp(w2)*trainData) + param.noiseSig*randn(param.nRecep,size(trainData,2));
MI = nonparanormal_info(resp');
H0 = 0;
for i0 = 1:param.nRecep
    H0 = H0 + HShannon_KDP_estimation(resp(i0,:),col);
end
H = -H0 + MI; 
disp(['the entropy for original optimal W is: ',num2str(H)])
% generate figureS2d
% this script calculates how the representation entropy changes when
% shuffled the optimal matrix
% this is used we we want to demonstrate that weak correlation among matrix
% elements contribute to the entropy

clear
clc

outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

%% define the param struct
nOdor = 100;
nRecp = 30;
spar = 2;
sig =2;
        
sparsity = spar;
        
param.nRecep = nRecp;
param.nOdor = nOdor;
param.lSig = sig;
param.lMu = 0;
param.nSamp = 2e4;
param.sparsity = true;
param.spRatio = sparsity/param.nOdor;
param.dType = 'lognorm';
param.eig = [];
param.corr = false;
param.regularize = false;
param.spType = 'absolute';          %type of sparisty, could be absolute or average
param.noiseSig = 1e-2;              % introduce small noise on the response
        
        
% parameter on W
% muW = -1;
% sigW = 2;
        
ShuffTimes = 200;    %randomly shuffling 100 times    
NW = 10;  % 10 times repeats
allH = zeros(ShuffTimes,NW);
allH0 = zeros(NW,1);

% load the example matrix
% dFolder = '/Users/shan/Dropbox/olfactionProject/data/GcmiDifSigData';
dFolder = '../data';

fName = 'N100_R30_S2_sig2_2018-05-18.mat';
load(fullfile(dFolder,fName));
w = reshape(allMat(:,10),[nRecp,nOdor]);


% initialize information estimator toolbox
mult = 1;
col = HShannon_KDP_initialization(mult);

% first, to see how the performance of optimal matrix
for i0 = 1:NW
    rng(100)
    eigVal = specifyEig(param.nOdor,param.eig);
    corrCoefMat = randCorrCoef('buildin',eigVal);
    trainData = genTrainData(param,corrCoefMat);
    
    % add some noise on the input o not
    if ~isempty(param.noiseSig)
       resp = reshape(exp(w),nRecp,nOdor)*trainData./(1+reshape(exp(w),nRecp,nOdor)*trainData)...
            + param.noiseSig*randn(param.nRecep,param.nSamp);
    else
       resp = reshape(exp(w),nRecp,nOdor)*trainData./(1+reshape(exp(w),nRecp,nOdor)*trainData);
    end
    
    % calculate the differential entropy for the representation
    MI = nonparanormal_info(resp');
        
        % using gcmi algorithm
    H0 = 0;
    for k0 = 1:param.nRecep
            H0 = H0 + HShannon_KDP_estimation(resp(k0,:),col);
    end
%         H = -H0 + MI - param.nRecep/2*log(2*pi*exp(1)*param.noiseSig^2);   %without noise
    H = -H0 + MI;   %without noise
    allH0(i0) = H;
end


% second, differential entropy after shuffling
for i0 = 1:ShuffTimes
    newW = reshape(w(randperm(nOdor*nRecp)),size(w,1),size(w,2));
    
    for j0 = 1:NW

%         trainData0 = zeros(param.nOdor,param.nSamp);
        eigVal = specifyEig(param.nOdor,param.eig);
        corrCoefMat = randCorrCoef('buildin',eigVal);
        trainData = genTrainData(param,corrCoefMat);
        
        % add some noise on the input o not
        if ~isempty(param.noiseSig)
            resp = reshape(exp(newW),nRecp,nOdor)*trainData./(1+reshape(exp(newW),nRecp,nOdor)*trainData)...
                + param.noiseSig*randn(param.nRecep,param.nSamp);
        else
            resp = reshape(exp(newW),nRecp,nOdor)*trainData./(1+reshape(exp(newW),nRecp,nOdor)*trainData);
        end
        
        % calculate the differential entropy for the representation
        MI = nonparanormal_info(resp');
        
        % using gcmi algorithm
        H0 = 0;
        for k0 = 1:param.nRecep
            H0 = H0 + HShannon_KDP_estimation(resp(k0,:),col);
        end
%         H = -H0 + MI - param.nRecep/2*log(2*pi*exp(1)*param.noiseSig^2);   %without noise
        H = -H0 + MI;   %without noise
        allH(i0,j0) = H;
    end
end

%% plot the figure
% looad data
load(fullfile(outFolder,'gcmi_diffEnt_OptW_randomShuffle_N100M30sig2sp2_03-Oct-2018_2.mat'))
% set graphics 
defaultGraphicsSetttings
%define some colors using brewermap
Bu = brewermap(11,'Blues');    % blues
Gr = brewermap(11,'Greys');    % greys

figure
hold on
histogram(-mean(allH,2),20,'Normalization','probability','FaceColor',Bu(9,:))
Yrange = get(gca,'Ylim');
plot(-[mean(allH0);mean(allH0)],Yrange','--')
legend('shuffle','original')
legend boxoff
box on
xlabel('differential entropy')
ylabel('probability')
prefix = ['gcmi_entrComp_OptW_Shuff','_N',num2str(nOdor),'M',num2str(nRecp),'sig',num2str(sig),....
    'sp',num2str(spar),'_',date,'_2'];
saveas(gcf,[outFolder,filesep,prefix,'.fig'])
print('-painters','-depsc',[outFolder,filesep,prefix,'.eps'])
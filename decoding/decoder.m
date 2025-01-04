%% training a neural network to perform decoding of OSN response
% this is a revision of Qianyi's original script
% version 0.4
% introduce noise in input
% last revised on 09/22/2018

function [] = decoder(nOdor,nRecp,spar,sig,noiseSig,varargin)

% parse input parameter, whether use random bumber seed
if nargin > 5
    rngSeed = varargin{1};
    rng(rngSeed)
else
    rng('shuffle')
end

if nargin > 6
    thd = varargin{2};
else
    thd = exp(-2*sig);% the detectable threshold 
end
% add search path to the working directory
addpath(genpath('/lustre1/tangc_pkuhpc/ssqin/olfactoryCoding/optiInterMatrix/myFunc'));
spAll = 0.05:0.05:0.95;

hiddenSize = 100;
% thd = exp(-2*sig);     

% random number generator seeds


% cycle through different number of hidden layers
Layers = 5;
tr_all = cell(length(spAll),Layers);
for j0 = 1:Layers
    for i0 = 1:length(spAll)
        sp = spAll(i0);
        %load('gcmi_N20_R9_S2_sig2_2018-03-20.mat')
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
        param.noiseSig = noiseSig;              % introduce small noise on the response

        
        % get the optimal mean and std of W
%         [muW,sigW,rhoW,fmin] = selectMeanSig(nOdor,spar);
        muW = -1;
        sigW = 2;
        w = normrnd(muW,sigW,[nRecp,nOdor]);
        w_zero = rand(nRecp,nOdor);
        w(w_zero>sp) = -inf;

%         trainData0 = zeros(param.nOdor,param.nSamp);
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

        %trainData0(trainData==0) = 0;
        %trainData0(trainData>0) = 1;
        %trainData = reshape(trainData(trainData>0),sparsity,param.nSamp);
        trainData = log(thd + trainData);

        % define the neural nework, with 100 units of hidden layer
        net = feedforwardnet(hiddenSize*ones(1,j0)); 
        %net.layers{end}.transferFcn = 'logsig';
        net.trainFcn = 'trainscg';
        net.trainParam.epochs = 1e4;
        net.performFcn = 'mse';  % mae, mse

        % train the network
        [net,tr] = train(net,resp,trainData);
        %net.performFcn = 'mse';
        %[net,tr] = train(net,resp,[trainData;trainData0]);
        tr_all{i0,j0}(:,1) = tr.best_perf;
        tr_all{i0,j0}(:,2) = tr.best_tperf;

%         outputs = net(resp);
        
%         plotFigure(outputs,trainData,tr,sp,param.noiseSig)
%         errors = gsubtract(trainData,outputs);
%         performance = perform(net,trainData,outputs);

        % visualize the traininig quality
%         tsOut = outputs(:,tr.testInd);
%         tsTarg = trainData(:,tr.testInd);
%         plotregression(tsTarg,tsOut,'Testing')
% 
%         temp1 = sort(tsOut,'descend');
%         temp2 = sort(tsTarg,'descend');
%         tsele = temp1(1:2,:);
%         trele = temp2(1:2,:);
        

%         bodyfatOutputs = net(bodyfatInputs);
%         trOut = bodyfatOutputs(tr.trainInd);
%         vOut = bodyfatOutputs(tr.valInd);
%         tsOut = bodyfatOutputs(tr.testInd);
%         trTarg = bodyfatTargets(tr.trainInd);
%         vTarg = bodyfatTargets(tr.valInd);
%         tsTarg = bodyfatTargets(tr.testInd);
%         plotregression(trTarg, trOut, 'Train', vTarg, vOut, 'Validation', tsTarg, tsOut, 'Testing')
    end
end
% save the final trainging results
save(strcat('nOdor',num2str(nOdor),'_nRecp',num2str(nRecp),'_sp',num2str(spar),'_hiddenSize',num2str(hiddenSize),'.mat'),'tr_all','spAll');
end

function plotFigure(outputs,trainData,tr,sp,noiseSig)
        
        % set the folder to save the figure
        saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

        % set the properties of graphics
        defaultGraphicsSetttings
        % first the regression figure, showing all the elements
        tsOut = outputs(:,tr.testInd);
        tsTarg = trainData(:,tr.testInd);
        plotregression(tsTarg,tsOut,'Testing')
        
        %second,only show the two selected, biggest
        temp1 = sort(tsOut,'descend');
        temp2 = sort(tsTarg,'descend');
        tsele = temp1(1:2,:);
        trele = temp2(1:2,:);
        
        figure
        hold on
        scatter(trele(:),tsele(:),10,'filled')
        xl = get(gca,'XLim');
%         yl = get(gca,'YLim');
        plot(xl,xl,'k--','LineWidth',1.5)
        hold off
        ylim(xl)
        
        xlabel('$\ln(1+x_i)$','Interpreter','latex')
        ylabel('$\ln(1+\hat{x}_i)$','Interpreter','latex')
        
        figNamePref = ['N20M9Sig2_recons_performance_spW',num2str(sp),'noiseSig',...
            num2str(noiseSig),'_',date];
        saveas(gcf,[saveFolder,filesep,figNamePref,'.fig'])
        print('-depsc',[saveFolder,filesep,figNamePref,'.eps'])
end
function [muW,sigW,rhoW,fmin] = selectMeanSig(N,sp)
% this function return the optimal parameter of parameterized W

% load the data
% dFoler = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';  % run local
dFoler = '../';  % run on cluster
dName= fullfile(dFoler,'gcmi_distri_summData_16-Jul-2018.mat');
load(dName)
allN = allN;
allSp = allSp;
inx1 = find(allN ==N);
inx2 = find(allSp ==sp);

muW = meanW(inx1,inx2);
sigW = meanSigW(inx1,inx2);
rhoW = meanRho(inx1,inx2);
fmin = meanFmin(inx1,inx2);
end
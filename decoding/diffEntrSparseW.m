% this script calculate how the representation entropy changes with the
% sparsity of W given the mu and sigma_w fixed 
% this is used we we want to demonstrate that maximum entropy coding
% facilitate downstream reconstruction

clear
clc

outFolder = '/home/shan/Documents/MATLAB/theoNeurosci/olfaction/decoding';

%% define the param struct
spAll = 0.05:0.05:1;
summH = zeros(length(spAll),2); %mean and std of entropy
for i0 = 1:length(spAll)
        sp = spAll(i0);
        %load('gcmi_N20_R9_S2_sig2_2018-03-20.mat')
        spar = 2;
        nRecp = 9;
        nOdor = 20;
        sig = 2;
        
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
        muW = -1;  %-1.78
        sigW = 2;
        
        % generate an instance of w
    NW = 20;  % 10 times repeats
    allH = zeros(NW,1);
    % initialize information estimator toolbox
    mult = 1;
    col = HShannon_KDP_initialization(mult);
        
    for j0 = 1:NW
        w = normrnd(muW,sigW,[nRecp,nOdor]);
        w_zero = rand(nRecp,nOdor);
        w(w_zero>sp) = -inf;

        trainData0 = zeros(param.nOdor,param.nSamp);
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
        allH(j0) = H;
    end
    summH(i0,:) = [mean(-allH),std(allH)];
%     mean(allH)
%     std(allH)
end

 %save the data
dName = ['entropyDistr_N',num2str(nOdor),'M',num2str(nRecp),'_sp',num2str(spar),...
   '_sig',num2str(sig),'_ns',num2str(param.noiseSig),'_',date,'.mat'];
save(fullfile(outFolder,dName),'spAll','summH')

% plot the sparsity dependent of entropy
% defaultGraphicsSetttings
Bu = brewermap(11,'Blues');    % blues

figure
hold on
errorbar(spAll',summH(:,1),summH(:,2),'o-','MarkerSize',10,'LineWidth',2,...
    'Color',Bu(10,:),'MarkerFaceColor',Bu(10,:))
[~,inx] = sort(summH(:,1),'descend');
ah = gca;
% Ymin = ah.YLim(1);
mqplot([spAll(inx(1));spAll(inx(1))],[ah.YLim(1),ah.YLim(2)],'k--','LineWidth',1)
hold off
lg = legend(['\mu =',num2str(muW) ,',\sigma_w =',num2str(sigW)],'Location','northwest');
set(lg,'FontSize',20)
xlabel('sparisty of W')
ylabel('differential entropy')
figNamePref = ['diffEntrSpW_fixMu',num2str(muW),'sig',num2str(sigW),'_N',...
    num2str(nRecp),'M',num2str(nOdor),'sp',num2str(spar),'_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])
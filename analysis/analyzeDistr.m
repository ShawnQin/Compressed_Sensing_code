% analyze the solution based on distribution: we simulataneous assume that
% the optimal matrix is approxmiated by a sparse lognormal distribution
% This can be compared with the full optimization method


close all
clear
clc

%% prepare the data
dataFolder = '/Users/shan/Dropbox/olfactionProject/data/gcmi/gcmi_distr';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
filesRaw = {allFile.name}';

% screen the files, only part are used. Since the folder may contain
% heterogeneous data set
str_marker = '2018-05-09';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

num = 40;  % number of repeat
h = 1;     % the Hill coefficient in the model

dataSumm = struct('fminMean',[],'fminStd',[],'muMean',[],'muStd',[],'spMean',[],...
    'spStd',[],'sigMean',[],'sigStd',[],'N',[],'allSig',[]);

% string flag to extract information
Spstr = '(?<= *_N)[\d.]+(?=_)';
sigStr = '(?<= *_sig)[\d.]+(?=_)';

allSp = [];
allSig = [];
for i0 = 1:length(files)
    allSp = [allSp,str2num(char(regexp(files{i0},Spstr,'match')))];
    allSig = [allSig, str2num(char(regexp(files{i0},sigStr,'match')))];
end

[val1,~,Rinx] = unique(allSp);
[val2,~,Sinx] = unique(allSig);
dataSumm.N = sort(val1);
dataSumm.allSig = sort(val2);

% add a matrix of cell array to save all the matrix
dataSumm.allW = cell(length(dataSumm.N),length(dataSumm.allSig));

dataSumm.fminMean = zeros(length(dataSumm.N),length(dataSumm.allSig));
dataSumm.fminStd = zeros(length(dataSumm.N),length(dataSumm.allSig));
for i0 = 1:length(files)
    temp = load(char(fullfile(dataFolder,filesep,files{i0})));
    if num ~= length(temp.allfmin)
        error('number of repeats in this simulation do not match!')
    else
        
    % exclude the outliner
    if Rinx(i0) == 4
       inx = 22;  %the data to be excluded
       temp.allfmin(inx) = nan;
       temp.allMat(:,inx) = nan;
    elseif Rinx(i0) == 5
       inx = 40;
       temp.allfmin(inx) = nan;
       temp.allMat(:,inx) = nan;     
    end
     
        sig = dataSumm.allSig(Sinx(i0));
        N = dataSumm.N(Rinx(i0));
        dataSumm.fminMean(Rinx(i0),Sinx(i0)) = nanmean(-temp.allfmin);
        dataSumm.fminStd(Rinx(i0),Sinx(i0)) = nanstd(temp.allfmin);
            
        dataSumm.muMean(Rinx(i0),Sinx(i0)) = nanmean(temp.allMat(1,:));
        dataSumm.muStd(Rinx(i0),Sinx(i0)) = nanstd(temp.allMat(1,:));
        
        dataSumm.sigMean(Rinx(i0),Sinx(i0)) = nanmean(temp.allMat(2,:));
        dataSumm.sigStd(Rinx(i0),Sinx(i0)) = nanstd(temp.allMat(2,:));
        
        dataSumm.spMean(Rinx(i0),Sinx(i0)) = nanmean(temp.allMat(3,:));
        dataSumm.spStd(Rinx(i0),Sinx(i0)) = nanstd(temp.allMat(3,:));
       
    end
 
end

%% plot figures
% first set the default figure parameters
defaultGraphicsSetttings

% define the colors
blues = brewermap(11,'Blues');
greys = brewermap(11,'Greys');

% plot N-dependent fmin
figure
errorbar(dataSumm.N',dataSumm.fminMean,dataSumm.fminStd,'o-','MarkerSize',12,'LineWidth',2)
xlabel('number of odorants')
ylabel('differential entropy')
legend('M = 20')
legend boxoff
figPref = ['Nx20_gcmiDistr_DiffEntr_sig_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


% plot N-dependence of mu
figure
errorbar(dataSumm.N',dataSumm.muMean,dataSumm.muStd,'o-','MarkerSize',12,'LineWidth',2)
ylim([-2,0])
xlabel('number of odorants')
ylabel('$\mu$','Interpreter','latex')
legend('M = 20')
legend boxoff
figPref = ['Nx20_gcmiDistr_mu_sig_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

%plot N-dependence of \sigma
figure
errorbar(dataSumm.N',dataSumm.sigMean,dataSumm.sigStd,'o-','MarkerSize',12,'LineWidth',2)
ylim([0,2])
xlabel('number of odorants')
ylabel('$\sigma_w$','Interpreter','latex')
legend('M = 20')
legend boxoff
figPref = ['Nx20_gcmiDistr_sigma_sig_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


%plot N-dependence of sparsity
figure
errorbar(dataSumm.N',dataSumm.spMean,dataSumm.spStd,'o-','MarkerSize',12,'LineWidth',2)
ylim([0,1])
xlabel('number of odorants')
ylabel('sparsity of W')
legend('M = 20')
legend boxoff
figPref = ['Nx20_gcmiDistr_sp_sig_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

%%  check the convergence of the method
% load the data
dfile = 'distr2_h2_N50_R10_S2_sig2_2018-05-12.mat';
load(fullfile(dataFolder,dfile));

%exclude bad data
goodInx = allfmin < 10;


% distribution of fmin
figure
histogram(-allfmin(goodInx),'Normalization','probability')
lg = legend('N = 20, M = 10, \sigma_c = 2, sp =2');
legend boxoff
set(lg,'FontSize',16)
xlabel('differentiatial entropy')
ylabel('probability')
figPref = ['50x10_gcmiDistr_fminHist','_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


%ditrbution mu
figure
histogram(allMat(1,goodInx),'Normalization','probability')
lg = legend('N = 20, M = 10, \sigma_c = 2, sp =2');
legend boxoff
set(lg,'FontSize',16)
xlabel('$\mu$','Interpreter','latex')
ylabel('probability')
figPref = ['50x10_gcmiDistr_muHist','_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


%distritbution of sigma_w
figure
histogram(allMat(2,goodInx),'Normalization','probability')
lg = legend('N = 20, M = 10, \sigma_c = 2, sp =2');
legend boxoff
set(lg,'FontSize',16)
xlabel('$\sigma_w$','Interpreter','latex')
ylabel('probability')
figPref = ['50x10_gcmiDistr_sigmaHist','_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


% distribution of sp
figure
histogram(allMat(3,goodInx),'Normalization','probability')
lg =  legend('N = 20, M = 10, \sigma_c = 2, sp =2');
legend boxoff
set(lg,'FontSize',16)
xlabel('sparsity of W')
ylabel('probability')
figPref = ['50x10_gcmiDistr_spHist','_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


% correlation among these parameters
set(groot,'DefaultFigurePaperUnits', 'inches', ...
           'DefaultFigureUnits', 'inches', ...
           'DefaultFigurePaperPosition', [0, 0, 10.1, 8.1], ...
           'DefaultFigurePaperSize', [10.1, 8.1], ...
           'DefaultFigurePosition', [0.1, 0.1, 10, 8]);
figure
hold on

% fmin-mu
C = corrcoef(allfmin(goodInx),allMat(1,goodInx));
subplot(2,3,1)
scatter(allMat(1,goodInx),allfmin(goodInx),8,blues(9,:),'filled')
legend(['\rho =',num2str(C(1,2))])
legend boxoff
xlabel('$\mu$','Interpreter','latex','FontSize',20)
ylabel('$f_{min}$','Interpreter','latex','FontSize',20)
set(gca,'FontSize',20)

% fmin - sigma_w
C = corrcoef(allMat(2,goodInx),allfmin(goodInx));
subplot(2,3,2)
scatter(allMat(2,goodInx),allfmin(goodInx),8,blues(9,:),'filled')
legend(['\rho =',num2str(C(1,2))])
legend boxoff
xlabel('$\sigma_w$','Interpreter','latex','FontSize',20)
ylabel('$f_{min}$','Interpreter','latex','FontSize',20)
set(gca,'FontSize',20)

% fmin - sp
C = corrcoef(allMat(3,goodInx),allfmin(goodInx));
subplot(2,3,3)
scatter(allMat(3,goodInx),allfmin(goodInx),8,blues(9,:),'filled')
legend(['\rho =',num2str(C(1,2))])
legend boxoff
xlabel('$sp$','Interpreter','latex','FontSize',20)
ylabel('$f_{min}$','Interpreter','latex','FontSize',20)
set(gca,'FontSize',20)


% mu-sigma
C = corrcoef(allMat(1,goodInx),allMat(2,goodInx));
subplot(2,3,4)
scatter(allMat(1,goodInx),allMat(2,goodInx),8,blues(9,:),'filled')
legend(['\rho =',num2str(C(1,2))])
legend boxoff
xlabel('$\mu$','Interpreter','latex','FontSize',20)
ylabel('$\sigma_w$','Interpreter','latex','FontSize',20)
set(gca,'FontSize',20)

% mu-sp
C = corrcoef(allMat(1,goodInx),allMat(3,goodInx));
subplot(2,3,5)
scatter(allMat(1,goodInx),allMat(3,goodInx),8,blues(9,:),'filled')
legend(['\rho =',num2str(C(1,2))])
legend boxoff
xlabel('$\mu$','Interpreter','latex','FontSize',20)
ylabel('$sp$','Interpreter','latex','FontSize',20)
set(gca,'FontSize',20)


% sigma_w - sp
C = corrcoef(allMat(2,goodInx),allMat(3,goodInx));
subplot(2,3,6)
scatter(allMat(2,goodInx),allMat(3,goodInx),8,blues(9,:),'filled')
legend(['\rho =',num2str(C(1,2))])
legend boxoff
xlabel('$\sigma_w$','Interpreter','latex','FontSize',20)
ylabel('$sp$','Interpreter','latex','FontSize',20)
set(gca,'FontSize',20)

hold off
figPref = ['50x10_gcmiDistr_paramCorr_','_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sparsity dependent effects
% unfortunatly, h = 2
%% prepare the data
%% prepare the data
dataFolder = '/Users/shan/Dropbox/olfactionProject/data/gcmi/gcmi_distr';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
filesRaw = {allFile.name}';

% screen the files, only part are used. Since the folder may contain
% heterogeneous data set
str_marker = '_N100_';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

num = 40;  % number of repeat
h = 2;     % the Hill coefficient in the model
sig = 2;   % standard deviation of input
N = 100;   % fixed the number of ligands

dataSumm = struct('fminMean',[],'fminStd',[],'muMean',[],'muStd',[],'spMean',[],...
    'spStd',[],'sigMean',[],'sigStd',[]);

% string flag to extract information
Spstr = '(?<= *_S)[\d.]+(?=_)';
RStr = '(?<= *_R)[\d.]+(?=_)';

allSp = [];
allR = [];
for i0 = 1:length(files)
    allSp = [allSp,str2num(char(regexp(files{i0},Spstr,'match')))];
    allR = [allR,str2num(char(regexp(files{i0},RStr,'match')))];
end
[val1,~,Rinx] = unique(allSp);
[val2,~,Sinx] = unique(allR);
dataSumm.allSp = sort(val1);
dataSumm.allR = sort(val2);

% add a matrix of cell array to save all the matrix

for i0 = 1:length(files)
    temp = load(char(fullfile(dataFolder,filesep,files{i0})));
    if num ~= length(temp.allfmin)
        error('number of repeats in this simulation do not match!')
    else
        
    % exclude the outliner
    inx = find(temp.allfmin > 1e2);
    temp.allfmin(inx) = nan;
     
        dataSumm.fminMean(Rinx(i0),Sinx(i0)) = nanmean(-temp.allfmin);
        dataSumm.fminStd(Rinx(i0),Sinx(i0)) = nanstd(temp.allfmin);
            
        dataSumm.muMean(Rinx(i0),Sinx(i0)) = nanmean(temp.allMat(1,:));
        dataSumm.muStd(Rinx(i0),Sinx(i0)) = nanstd(temp.allMat(1,:));
        
        dataSumm.sigMean(Rinx(i0),Sinx(i0)) = nanmean(temp.allMat(2,:));
        dataSumm.sigStd(Rinx(i0),Sinx(i0)) = nanstd(temp.allMat(2,:));
        
        dataSumm.spMean(Rinx(i0),Sinx(i0)) = nanmean(temp.allMat(3,:));
        dataSumm.spStd(Rinx(i0),Sinx(i0)) = nanstd(temp.allMat(3,:));
       
    end
end

dataName = fullfile(saveFolder,'gcmi_distr_h2_diffSp.mat');
save(dataName,'dataSumm')

%% plot the figures

% plot the differential entropy as a function of the sparsity
figure
errorbar(dataSumm.allSp'*ones(1,2),dataSumm.fminMean,dataSumm.fminStd,'o-',...
    'MarkerSize',12,'LineWidth',2)
legend('R = 10','R=20')
legend boxoff
xlabel('number of odorants')
ylabel('differential entropy')
figPref = ['gcmi_distr_diffEntropy_sp_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

% plot the mean Wij as a function of the sparsity
figure
errorbar(dataSumm.allSp'*ones(1,2),dataSumm.muMean,dataSumm.muStd,'o-',...
    'MarkerSize',12,'LineWidth',2)
legend('R = 10','R=20')
legend boxoff
xlabel('number of odorants')
ylabel('average of $W_{ij}$','Interpreter','latex')
figPref = ['gcmi_distr_mu_sp_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

% plot the sparsity of W as a function of the sparsity
figure
errorbar(dataSumm.allSp'*ones(1,2),dataSumm.spMean,dataSumm.spStd,'o-',...
    'MarkerSize',12,'LineWidth',2)
legend('R = 10','R=20')
legend boxoff
xlabel('number of odorants')
ylabel('sparity of $W_{ij}$','Interpreter','latex')
figPref = ['gcmi_distr_sparsity_sp_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

% plot the std of W as a function of the sparsity
figure
errorbar(dataSumm.allSp'*ones(1,2),dataSumm.sigMean,dataSumm.sigStd,'o-',...
    'MarkerSize',12,'LineWidth',2)
legend('R = 10','R=20')
legend boxoff
xlabel('number of odorants')
ylabel('$\sigma_{w}$','Interpreter','latex')
figPref = ['gcmi_distr_sigW_sp_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])



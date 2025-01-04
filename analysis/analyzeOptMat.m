% analyze and plot the optimal interaction matrix

close all
clear
clc

% dataFolder = '/home/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/spLogN50_R20_S2_sig3';
dataFolder = '/home/shan/Documents/MATLAB/theoNeurosci/olfaction/spLogN100_R16_S2_sig2_2017-12-26';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
files = {allFile.name}';

%merge data
numOdor = 100;
numRecp = 16;

sparsity =2; %average number of ligand in mixture
numSamp = 1e4; 
mu = 0;
sig = 2; 

%%
allMat = zeros(numOdor*numRecp,length(files));
allfmin = zeros(length(files),1);
for i0 = 1:length(files)
    temp = load(char(fullfile(dataFolder,filesep,files{i0})));
    allMat(:,i0)  = temp.wmin;
    allfmin(i0) = temp.fmin;
end

dateInfo = datestr(now,'yyyy-mm-dd');
% save(['N',num2str(numOdor),'_R',num2str(numRecp),'_sp',num2str(sparsity),...
%     '_sig',num2str(sig),'_',dateInfo,'.mat'],'allMat','allfmin')
save(['N',num2str(numOdor),'_R',num2str(numRecp),...
    '_sig',num2str(sig),'_',dateInfo,'.mat'],'allMat','allfmin')
%%
%plot the overall elements of optimal matrix
figure(1)
ROW = 6;
COL = 9;
for j0 = 1:ROW
    for k0 = 1:COL
        subplot(ROW,COL,(j0-1)*COL+k0)
        histogram(allMat(:,(j0-1)*COL+k0),20)
        set(gca,'XLim',[-20 5])
    end
end

for j0 = 1:size(allMat,2)
    subplot(ROW,COL,j0)
     histogram(allMat(:,j0),30)
     set(gca,'XLim',[-20 20])
end


%overall distribution
figure(2)
histogram(allMat(abs(allMat)<10),50)

% histogram of output
mult = 1;
col1 = HShannon_KDP_initialization(mult);
col2 = HShannon_vME_initialization(mult);
wmin = reshape(allMat(:,1),[numRecp,numOdor]);
trainData = zeros(numOdor,numSamp);
inx = datasample(1:numOdor*numSamp,numSamp*sparsity,'Replace',false);
trainData(inx) = exp(normrnd(mu,sig,[numSamp*sparsity,1]));  
h1 = exp(wmin) * trainData./(1+ exp(wmin) * trainData);
HShannon_KDP_estimation(h1,col1)
% HShannon_vME_estimation(h1,col2)
figure
for i0 = 1:numRecp
    subplot(1,numRecp,i0)
    histogram(h1(i0,:),50)
end

%effective dimension
C = cov(h1');
effDim = trace(C)^2/trace(C*C);

%% analysis synthesis data
% load the data

dName = 'N50_R17_sp3_sig2_2017-12-02.mat';
load(dName)

N=50;
M=17;
S=3;
sig=2;
titleName=['N=',num2str(N),',R=',num2str(M),',S=',num2str(S),',sig=',num2str(sig)];
map = brewermap(9,'RdBu'); 

% first plot the over all histogram
figure
% histogram(allMat(allMat>-20),50)
histogram(allMat,100,'Normalization','probability')
title(titleName)
xlim([-80 15])
xlabel('log(Wij)','FontSize',24)
ylabel('probability','FontSize',24)
set(gca,'LineWidth',1.5,'FontSize',20)


figure
fh = histfit(allMat(allMat>-80),50,'kernel');
fh(1).FaceColor = map(8,:);
fh(1).EdgeColor = [0 0 0];
xlim([-80 15])
title(titleName)
xlabel('log(Wij)','FontSize',24)
ylabel('Count','FontSize',24)
set(gca,'LineWidth',1.5,'FontSize',20)


% plot individual
inxSel = 1:6;
figure
title(titleName)
for i0 = 1:6
    subplot(2,3,i0)
    temp = allMat(:,inxSel(i0));
    histogram(temp(temp > -20),'Normalization','probability')
    xlim([-20,15])
    xlabel('log(Wij)','FontSize',18)
    ylabel('probability','FontSize',18)
    set(gca,'LineWidth',1,'FontSize',16)
end

%% this part try to visualize diffent optimal matrix
dataFolder = '/Users/shan/Dropbox/olfactionProject/data/optMatAnal';

% all file related to R = 20, Sigma=1.5 ,S=2
allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
allFile = {allFile.name}.';  %get all the name of mat files

%get indices of files that match regular expression
% define the function
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(allFile,str,'once'));
pattern = '_R20_sp2_sig1.5_';
allSel =  allFile(FIND(pattern));

%sort N
% SORT = @(str,cellArray) regexp(cellArray,str,'match');
str1 = '^N(\d+)_*';
% ext = regexp(temp, str1, 'tokens', 'once');
% extractN = @(s) cellfun(@(s) ~isempty(s), regexp(s, str1, 'tokens', 'once'));
% allN = allSel(extractN(str1));
allN = nan(length(allSel),1);
Wselect = cell(length(allSel),1);
allData = cell(length(allSel),1);
for i0 = 1:length(allSel)
    dataFile = fullfile(dataFolder,filesep,allSel{i0});
    load(dataFile)
    allData{i0} = allMat;
    allN(i0) = str2num(char(regexp(allSel{i0}, str1, 'tokens', 'once')));
    col = size(allMat,2);
    Wselect{i0} = reshape(allMat(:,randi(col)),[20,allN(i0)]);
%     figure
%     set(gcf,'Units','inches','Position',[4,4,2,3.8])
%     imagesc(Wselect{i0},[-5,5])
%     colormap(jet)
%     title(['N=',num2str(allN(i0))])
end

% plot all the histogram
% [ha, pos] = tight_subplot(2,2,[.05 .1],[.2 .05],[.15 .05]);
figure
[val,inx] = sort(allN);
for i0 = 1:length(inx)
    dataFile = fullfile(dataFolder,filesep,allSel{inx(i0)});
    load(dataFile)
%     axes(ha(i0))
%     hold on
    subplot(3,5,i0)
    histogram(allMat(allMat>-20),40)
%     legend(['N=',num2str(allN(i0))])
%     legend boxoff
    title(['N=',num2str(val(i0))],'FontSize',14)
    xlim([-20,10])
    set(gca,'LineWidth',1,'FontSize',14)
    if ismember(i0, [1,6,11])
        ylabel('Count','FontSize',16)
    else
%         yticks([])
    end
    
    if i0 > 11
        xlabel('log(W)','FontSize',16)
    else
%         xlabel([])
    end
%     hold off
end

selectN = [3,4,5,8,10,12,15,20,30,50,80];


% plot heatmap of the matrix, use tight plot

%% this part analyze the rank correlation
% using Spearman rank correlation to see if there is a transition in the
% structure of optimal matrix when increasing N
[val,inx] = sort(allN);
allOdorCorr = cell(length(allSel),1);
allShufCorr = cell(length(allSel),1);
numRecp = 20;   %number of receptors
sig = 1.5;
LB = -4*sig; %lower bound, smaller than which is considered no response
UB = 4*sig;  % uper bound, larger than which is considered to be satureated
for i0 = 1:length(inx)
    dataFile = fullfile(dataFolder,filesep,allSel{inx(i0)});
    load(dataFile)
    % loop through different trials, or different matrix
    allOdorCorr{i0} = nan(allN(inx(i0))*(allN(inx(i0))-1)/2,size(allMat,2));
    allShufCorr{i0} = nan(allN(inx(i0))*(allN(inx(i0))-1)/2,size(allMat,2));
    for l0 = 1:size(allMat,2)
        temp = reshape(allMat(:,l0),[numRecp,allN(inx(i0))]);
        temp(temp < LB) = LB;
        temp(temp > UB) = UB;
        
        % randomly shuffled matrix
        NUM = numRecp*allN(inx(i0));
        temp2 = reshape(temp(randperm(NUM)),[numRecp,allN(inx(i0))]);
        
         % loop through different pairs of columns
        corrRecord = [];
        for j0 = 1:(allN(inx(i0))-1)
            for k0 = (j0+1):allN(inx(i0))
                nonzeroInx = any(temp(:,[j0,k0]) ~= LB,2); % not all are zero
%                 saturateInx  = any(temp(:,[j0,k0]) == UB);
                selectColEle = temp(nonzeroInx,[j0,k0]);
                corrRecord = [corrRecord,corr(selectColEle(:,1),selectColEle(:,2),'type','Spearman')];
                
                nonzeroInx2 = any(temp2(:,[j0,k0]) ~= LB,2); % not all are zero
                selectColEle2 = temp(nonzeroInx2,[j0,k0]);
                corrRecord2 = [corrRecord2,corr(selectColEle2(:,1),selectColEle2(:,2),'type','Spearman')];
            end
        end
        allOdorCorr{i0}(:,l0) = corrRecord';
        allShufCorr{i0}(:,l0) = corrRecord2';
    end  
end


% plot the rank correlation distribution
figure
% set(gca,'Units','inches','Position',[0,0,8,6])
hold on

% set the color theme
allColor = brewermap(15,'Spectral');
corrMean = zeros(length(allN),1);
corrStd = zeros(length(allN),1);
for i0 = 1:length(allN)
    corrMean(i0) = mean(allOdorCorr{i0}(:));
    corrStd(i0) = std(allOdorCorr{i0}(:));
    [fi,xi] = ksdensity(allOdorCorr{i0}(:));
    plot(xi,fi,'LineWidth',2,'Color',allColor(i0,:),'DisplayName',['N=',num2str(allN(inx(i0)))])
end
hold off
set(gca,'FontSize',24,'LineWidth',1)
xlabel('Spearman rank correlation','FontSize',30)
ylabel('Probability','FontSize',30)


figure
errorbar(allN(inx),corrMean,corrStd)
set(gca,'LineWidth',1,'FontSize',24)
xlabel('Number of ligands','FontSize',28)
ylabel('Mean rank correlation','FontSize',28)

% for the randomly shuffled counterpart
figure
hold on
corrMeanShuf = zeros(length(allN),1);
corrStdShuf = zeros(length(allN),1);
for i0 = 1:length(allN)
    corrMeanShuf(i0) = mean(allShufCorr{i0}(:));
    corrStdShuf(i0) = std(allShufCorr{i0}(:));
    [fi,xi] = ksdensity(allShufCorr{i0}(:));
    plot(xi,fi,'LineWidth',2,'Color',allColor(i0,:),'DisplayName',['N=',num2str(allN(inx(i0)))])
end

figure
errorbar(allN(inx),corrMeanShuf,corrStdShuf)
set(gca,'LineWidth',1,'FontSize',24)
xlabel('Number of ligands','FontSize',28)
ylabel('Mean rank correlation','FontSize',28)


%% This part analyzes the 2x2 simulation data

% dataFolder = './';

allSig = 0.5:0.5:4;

sumDiag = zeros(length(allSig),2);    % mean and std
sumOffDiag = zeros(length(allSig),2); % mean and std
allEntropy = zeros(length(allSig),2); % mean and std target function
for i0 = 1:length(allSig)
    fileName = fullfile(['N2_R2_sp2_sig',num2str(allSig(i0)),'_2017-12-14.mat']);
    load(fileName)
    temp = sort(allMat);
    sumDiag(i0,1) = mean([temp(1,:),temp(2,:)]/allSig(i0));
    sumDiag(i0,2) = std([temp(1,:),temp(2,:)]/allSig(i0));
    sumOffDiag(i0,1) = mean([temp(3,:),temp(4,:)]/allSig(i0));
    sumOffDiag(i0,2) = std([temp(3,:),temp(4,:)]/allSig(i0));
    allEntropy(i0,:) = [mean(allfmin),mean(allfmin)];
end

% plot the sigma dependent off diagnal and diagnoal elements
allColor = brewermap(15,'Spectral');  %define the colors
figure;
errorbar(allSig',sumDiag(:,1),sumDiag(:,2),'o-','MarkerSize',12,'LineWidth',2)
hold on
errorbar(allSig',sumOffDiag(:,1),sumOffDiag(:,2),'o-','MarkerSize',12,'LineWidth',2)
hold off
% set(fh,{'DisplayName'},{'diagonal';'off diagonal'})
% legend show
% legend('show')
xlim([0 5])
ylim([-30 1])
xlabel('\sigma','FontSize',30)
ylabel('value/\sigma','FontSize',30)
set(gca,'LineWidth',1,'FontSize',24)


% plot the fitness function
figure
errorbar(allSig',allEntropy(:,1),allEntropy(:,2),'o-','MarkerSize',15,'LineWidth',3)
xlim([0 5])
xlabel('\sigma','FontSize',30)
ylabel('cost function','FontSize',30)
set(gca,'LineWidth',1,'FontSize',24)

%% this part analyze 2 by M data
% we compare the target function, the rank correlation and the average
% elements values

%load the data
% allR = [2 3 4 5 10 15 20 50 100];
allR = [1:1:20,30,40,50,100,200];

sigma = 1;
sumFmin = zeros(length(allR),2);   % mean and std of fmin
sumActive = zeros(length(allR),2); % mean and std of active elements
sumSparsity = zeros(length(allR),2); %mean and average active elements
% sumRankCorr = zeros(sum(allR>=10),2);  %rank correlation
sumRankCorr = [];  %rank correlation

for i0 = 1:length(allR)
    fileName = fullfile(['N2_R',num2str(allR(i0)),'twoByM_sig1_2018-01-27']);
    load(fileName)
    sumFmin(i0,:) = [-mean(allfmin) + log(2*pi*exp(1)*sigma^2),std(allfmin)];
%     temp = allMat;
%     temp(allMat<4*sigma) = nan;  %set nonresponse elements as nan
    sumActive(i0,:) = [mean(allMat(allMat > -4*sigma)),std(allMat(allMat > -4*sigma))];
%     temp2 = sum(allMat>-4*sigma);
    sumSparsity(i0,:) = [mean(sum(allMat>-4*sigma)/allR(i0))/2,std(sum(allMat>-4*sigma)/allR(i0)/2)];
    nrec = size(allMat,1)/2;
    temp = zeros(length(allfmin),1);
    if allR(i0) >=5
    for j0 = 1:length(allfmin)
        temp(j0) = corr(allMat(1:nrec,j0),allMat((1+nrec):end,j0),'type','Spearman');
    end
    sumRankCorr = [sumRankCorr;[mean(temp),std(temp)]];
    end
end

% load another data set, different sigma
sumFmin2 = zeros(length(allR),2);   % mean and std of fmin
sumActive2 = zeros(length(allR),2); % mean and std of active elements
sumSparsity2 = zeros(length(allR),2); %mean and average active elements
sumRankCorr2 = [];  %rank correlation
sigma = 2;
for i0 = 1:length(allR)
    fileName = fullfile(['N2_R',num2str(allR(i0)),'_sig2_2017-12-16.mat']);
    load(fileName)
    sumFmin2(i0,:) = [-mean(allfmin) + log(2*pi*exp(1)*sigma^2),std(allfmin)];
%     temp = allMat;
%     temp(allMat<4*sigma) = nan;  %set nonresponse elements as nan
    sumActive2(i0,:) = [mean(allMat(allMat > -4*sigma)),std(allMat(allMat > -4*sigma))];
    temp2 = sum(allMat>-4*sigma);
    sumSparsity2(i0,:) = [mean(sum(allMat>-4*sigma)/allR(i0))/2,std(sum(allMat>-4*sigma)/allR(i0)/2)];
    nrec = size(allMat,1)/2;
    temp = zeros(length(allfmin),1);
    if allR(i0) >=5
    for j0 = 1:length(allfmin)
        temp(j0) = corr(allMat(1:nrec,j0),allMat((1+nrec):end,j0),'type','Spearman');
    end
    sumRankCorr2 = [sumRankCorr2;[mean(temp),std(temp)]];
    end
end

% plot the figures
rdColor = brewermap(11,'RdBu');
figure
fh1a = errorbar(allR',sumFmin2(:,1),sumFmin2(:,2),'o-','MarkerSize',12,'LineWidth',2,'Color',rdColor(2,:));
hold on
fh1b = errorbar(allR',sumFmin(:,1),sumFmin(:,2),'o-','MarkerSize',12,'LineWidth',2,'Color',rdColor(9,:));
legend([fh1a fh1b],{'\sigma=2','\sigma=4'})
hold off
% legend(['\sigma =',num2str(sigma)])
xlabel('number of receptors','FontSize',28)
ylabel('diffeential entropy','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
saveas(gcf,'2xM_fmin.fig')

% plot rank correlation
figure
temp2 = allR(allR>=5);
errorbar(temp2',sumRankCorr2(:,1),sumRankCorr2(:,2),'o-','MarkerSize',12,'LineWidth',2,'Color',rdColor(2,:));
hold on
f2b = errorbar(temp2',sumRankCorr(:,1),sumRankCorr(:,2),'o-','MarkerSize',12,'LineWidth',2,'Color',rdColor(9,:));
hold off
legend([f2a f2b],{'\sigma=2','\sigma=4'})

ylim([-1.2,0])
xlabel('number of receptors','FontSize',28)
ylabel('rank correlation','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
saveas(gcf,'2xM_rankCorr.fig')

% plot mean value of active elements
figure
% yyaxis left
f3a = errorbar(allR,sumActive2(:,1),sumActive2(:,2),'o-','MarkerSize',12,'LineWidth',2,'Color',rdColor(2,:));
hold on
f3b = errorbar(allR,sumActive(:,1),sumActive(:,2),'o-','MarkerSize',12,'LineWidth',2,'Color',rdColor(9,:));
hold off
legend([f3a f3b],{'\sigma=2','\sigma=4'})

ylabel('mean value','FontSize',28)
% yyaxis right
xlabel('number of receptors','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
saveas(gcf,'2xM_meanVal.fig')

% plot sparsity
figure
f4a = errorbar(allR,sumSparsity2(:,1),sumSparsity2(:,2),'o-','MarkerSize',12,'LineWidth',2,'Color',rdColor(2,:));
hold on
f4b = errorbar(allR,sumSparsity(:,1),sumSparsity(:,2),'o-','MarkerSize',12,'LineWidth',2,'Color',rdColor(9,:));
legend([f4a f4b],{'\sigma=2','\sigma=4'})
hold off
xlabel('number of receptors','FontSize',28)
ylabel('sparsity','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
saveas(gcf,'2xM_sparsity.fig')


%% this part extract all the one by M data into concrete form
dataFolder = '/home/shan/Documents/MATLAB/theoNeurosci/olfaction/twoByM_sig1_18-01-27';
sig = 1;
allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
files = {allFile.name}';


allM = [1:20,30,40,50,100,200];
summFmin = [];
summWmin = cell(length(allM),1);

for i0 = 1:length(files)
    load(char(fullfile(dataFolder,filesep,files{i0})));
    summFmin = [summFmin,allFmin];
    for j0 = 1:length(allM)
        summWmin{j0} = [summWmin{j0},allWmin{j0}];
    end
end

dateInfo = datestr(now,'yyyy-mm-dd');
% save(['oneByM_sig',num2str(sig),'_',dateInfo,'.mat'],'summFmin','summWmin')
save(['twoByM_sig',num2str(sig),'_',dateInfo,'.mat'],'summFmin','summWmin')

mf = mean(summFmin,2);
stdf = std(summFmin,0,2);
errorbar(allM',mf,stdf)

histogram(summFmin(15,:))

% figure
% for i0 = 1:9
%     subplot(3,3,i0)
%     histogram(resp(i0,:))
%     xlim([0 1])
% end
% 
% figure
% count = 1;
% for i0 = 1:8
%     for j0 = i0:9
%         subplot(6,6,count)
%         scatter(resp(i0,randperm(1e4,2e3)),resp(j0,randperm(1e4,2e3)),10,'.')
%         count = count + 1;
%     end
% end
       

%% Plot data for one odor and many recptors
clear

% data folder
% dFolder = '/Users/shan/Dropbox/olfactionProject/data/oneByM01142018';
dFolder = '/Users/shan/Dropbox/olfactionProject/data/oneByM01142018';
% dFolder = './oneByM01142018';

% output folder, save the figures
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
allFile = dir(dFolder);
NUM = length(allFile) - 2;  %number of different data, for different sigma
allSig = zeros(1,length(NUM));
allM = [1:20,30,40,50,100,200];  % all different number of receptors
allFmin = zeros(length(allM),NUM); % store the average number of target function
allFminStd = zeros(length(allM),NUM);
allWmean = zeros(length(allM),NUM); %store the mean of Wij
allWstd = zeros(length(allM),NUM);  %store the std of Wij

pt = '(?<=_sig)[\d\.]+(?=_)';  % pattern to extract the value sigma
for i0 = 1:NUM
    fileName = allFile(2+i0).name;
    allSig(i0) = str2double(char(regexp(fileName,pt,'match')));
    load([dFolder,filesep,fileName]);
    allFmin(:,i0) = mean(summFmin(1:length(allM),:),2);
    allFminStd(:,i0) =  std(summFmin(1:length(allM),:),[],2);
    for j0 = 1:length(allM)
        temp = summWmin{j0}(:);
%         temp(temp>4)=4;   %out linear is set to the upper limit
%         temp(temp<-4)=-4; %out linear is set to the lower limit
        allWmean(j0,i0) = mean(temp(:));
        allWstd(j0,i0) = std(temp(abs(temp) <=4));
    end
end

% theoretical std of wij
h = 1;
orderSig = sort(allSig);
theoStd = sqrt(max(orderSig.^2 - 1/h^2,0));  % when sigma_c is smaller than 1

%plot the and compare with theory
myColor = brewermap(11,'Blues');

% differential entropy
figure
hold on
largeMeanStd = zeros(2,NUM); %store the mean and std of last data, very large M
for i0 = 1:NUM
    errorbar(allM',allFmin(:,i0),allFminStd(:,i0),'LineWidth',2)
    largeMeanStd(1,i0) = allFmin(end,i0);
    largeMeanStd(2,i0) = allFminStd(end,i0);
end
hold off
xlabel('M','FontSize',28)
ylabel('differential entropy','Interpreter','latex','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1)
figNamePref = 'oneByMdiffEntropySimu';
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.pdf'])

%plot differential entropy as a function of sigma_c
[~,Ix] = sort(allSig);
figure 
errorbar(orderSig,largeMeanStd(1,Ix),largeMeanStd(2,Ix),'LineWidth',2)
xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('differential entropy','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1,'XScale','log')
figNamePref = 'oneByMdiffEntropyLargeM';
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.pdf'])

% compare the wij with theory
figure
hold on
allStdOrder = sort(allWstd(25,:).*allSig);
plot(theoStd,allStdOrder,'o','MarkerEdgeColor',myColor(9,:),...
    'MarkerFaceColor',myColor(9,:),'LineWidth',2,'MarkerSize',15)
plot([0,10],[0,10],'k--','LineWidth',2)
hold off
xlim([0 11])
ylim([0 11])
xlabel('$\sqrt{\sigma_c^2 - 1/h^2}$','Interpreter','latex','FontSize',28)
ylabel('$\sigma(W)$','Interpreter','latex','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1)
figNamePref = 'oneByMCompSimuTheory';
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.pdf'])


S = 0:0.1:10;
theorSig = sqrt(max(S.^2 - 1/h^2,0));
figure
hold on
plot(S,theorSig,'k-','LineWidth',2)
plot(sort(allSig),allStdOrder,'o','MarkerEdgeColor',myColor(9,:),...
    'MarkerFaceColor',myColor(9,:),'MarkerSize',8)
xlim([-0.2 11])
ylim([-0.2,11])
xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('$\sigma(W)$','Interpreter','latex','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1)
figNamePref = 'oneByMsigWSimu';
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.pdf'])

% ===========================================================
% plot the position of optimal W along the cdf of input data
% ===========================================================
grayColor = brewermap(15,'RdGy');
sigSlect = 2;
cdfPosi = cell(length(allM),1);
figure
hold on
for i0 = 1:length(allM)
    temp = sort(summWmin{i0}(:,10),1);
    selectW = temp(abs(temp) < 5);
    cdfPosi{i0} = normcdf(selectW);
    plot(allM(i0),cdfPosi{i0},'ko','MarkerSize',6,'MarkerEdgeColor','k',...
        'MarkerFaceColor','k','LineWidth',1)
end
% hold off
xlabel('M','FontSize',28)
ylabel('cdf(W)','FontSize',28,'FontName','Helvetica')
set(gca,'FontSize',24,'LineWidth',1)
figNamePref = 'cdfIllus';
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.pdf'])

% figure
% hold on
allLine = cell(10,1);
% X = 0:0.2:20;
for i0  = 1:10
    allLine{i0} = [];
    for j0 = (2*i0-1):20  %only plot the first 20
        allLine{i0} = [allLine{i0},cdfPosi{j0}(i0)];
    end 
    X  = allM(2*i0-1):0.5:20;
    yy = spline(allM((2*i0-1):20),allLine{i0},X);
    plot(X,yy,'Color',grayColor(12,:),'LineWidth',1.5)
end

% plot the position on the coordinate of input
stimuPosi = cell(length(allM),1);
figure
hold on
for i0 = 1:length(allM)
    temp = sort(summWmin{i0}(:,10),1);
    stimuPosi{i0} = temp(abs(temp) < 5);
    plot(allM(i0),stimuPosi{i0},'ko','MarkerSize',5,'MarkerEdgeColor','k',...
        'MarkerFaceColor','k','LineWidth',1)
end
% hold off
xlabel('M','FontSize',28)
ylabel('W','FontSize',28,'FontName','Helvetica')
set(gca,'FontSize',24,'LineWidth',1)


% =======================
% the average position interval, to see if they are uniform in the large M
% limit
% ========================
edge = [0,1];  %edge position, since it is a cdf
meanInterval = zeros(length(allM),1);
stdInteval = zeros(length(allM),1);
for i0 = 1:length(allM)  %M = 1 only have one point
    temp = [edge(1);cdfPosi{i0};edge(2)];
    meanInterval(i0) = mean(diff(temp));
    stdInteval(i0) = std(diff(temp));
end

figure
X = 1./(allM+1);
errorbar(1./(allM+1),meanInterval,stdInteval,'o','MarkerSize',8,'LineWidth',2)
hold on
plot([0,0.5],[0,0.5],'k--','LineWidth',1.5)
hold off
xlabel('$1/(1+M)$','Interpreter','latex','FontSize',28)
ylabel('distance between adjacent Wi','Interpreter','latex','FontSize',28)
set(gca,'LineWidth',1,'FontSize',24)
figNamePref = 'oneByWinterval';
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.pdf'])


% ===========================================
% illustrate of pdf and cdf, using tight plot
% ===========================================

pdfGauss = chebfun('normpdf(x)',[-5,5]);
cdfGauss = chebfun('normcdf(x)',[-5,5]);

figure
set(gcf,'Units','inches','Position',[0,0,5,5],...
'PaperPositionMode','auto');
[ha, pos] = tight_subplot(2,1,[.02 .1],[.25 .05],[.2 .03]);

axes(ha(1))
hold on
plot(pdfGauss,'LineWidth',2)
ylabel('pdf','FontSize',24)
xticks = [];
set(gca,'LineWidth',1,'FontSize',20)

allX = norminv(0.2:0.2:0.8);
for i0 = 1:4
    plot([allX(i0),allX(i0)],[0,normpdf(0)],'--','Color',grayColor(12,:),'LineWidth',1)
end

axes(ha(2))
hold on
plot(cdfGauss,'LineWidth',2)
ylabel('cdf','FontSize',24)
set(gca,'XTick',-5:5:5,'LineWidth',1,'FontSize',20)

% add horizontal dashed lines
for i0 = 1:4
    plot([-5,5],[0.2*i0,0.2*i0],'--','Color',grayColor(12,:),'LineWidth',1)
end

sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
fileNameEps = [sFolder,filesep,'normpdf_cdf.eps'];
fileNameFig = [sFolder,filesep,'normpdf_cdf.fig'];
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)


%% The gradient-based method
histogram(log(allW(:)),'Normalization','probability')
xlabel('$\ln(W)$','Interpreter','latex','FontSize',28)
ylabel('probability','FontSize',28)
xlim([-20,10])
set(gca,'XTick',-20:5:10,'FontSize',24,'LineWidth',1.5)


finalVal = zeros(20,1);
for i0 = 1:20
    finalVal(i0) = totalCost{i0}(end);
end
meanF = mean(finalVal);
stdF = std(finalVal);

histogram(finalVal,'Normalization','probability')
xlabel('$\frac{1}{2}\langle \ln \det(\chi^T\chi)\rangle$','Interpreter','latex','FontSize',28)
ylabel('probability','FontSize',28)
set(gca,'XTick',-1.75:0.1:-1.5,'FontSize',24,'LineWidth',1.5)


%% This part is for the power law distribution of input
% dataFolder = '/Users/shan/Dropbox/olfactionProject/data/powerLaw';
dataFolder = '/Users/shan/Documents/olfactionData/powerLaw';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
files = {allFile.name}';


s1 = '(?<= *N)[\d.]+(?=_)';
s2 = '(?<= *_R)[\d]+(?=_)';
s3 = '(?<= *_S)[\d]+(?=_)';
s5 = '(?<= *sig)[\d.]+(?=_)';

numOdor = zeros(length(files),1);
numRecp = zeros(length(files),1);
spInput = zeros(length(files),1);
sig = zeros(length(files),1);
spW = zeros(length(files),1);

% figure
for i0 = 1:length(files)
    load(fullfile(dataFolder,files{i0}))
    numOdor(i0) = str2num(char(regexp(files{i0},s1,'match')));
    numRecp(i0) = str2num(char(regexp(files{i0},s2,'match')));
    spInput(i0) = str2num(char(regexp(files{i0},s3,'match')));
    sig(i0) = str2num(char(regexp(files{i0},s5,'match')));
    
    t = allMat(:);
%     subplot(3,3,i0)
%     
% %     histogram(t,'Normalization','probability')
%     histogram(t(abs(t)<20),'Normalization','probability')
% 
%     lg = ['N=',num2str(numOdor(i0)),', R = ',num2str(numRecp(i0)),', sp = ',num2str(spInput(i0)),...
%         ', \sigma = ',num2str(sig(i0))];
%     lgh = legend(lg,'Interpreter','latex');
%     set(lgh,'FontSize',14)
%     legend boxoff
%     xlim([-20,5])
%     set(gca,'FontSize',16)
%     xlabel('$\ln(w)$','Interpreter','latex','FontSize',20)
%     ylabel('probability','FontSize',20)
    spW(i0) = sum(t>-15)/length(t);
end

% group the data
uniqOdor = sort(uniq(numOdor));
uniqSp = sort(uniq(spInput));
sparsity = cell(length(uniqOdor),1);
spOut = cell(length(uniqOdor),1);
for i0 = 1:length(uniqOdor)
    inx = numOdor == uniqOdor(i0);
    sparsity{i0} = spInput(inx);
    spOut{i0} = spW(inx);
end

%plot the output sparsity and input sparisty
defaultGraphicsSetttings
figure
hold on
for i0 = 1:length(uniqOdor)
    plot(sparsity{i0},spOut{i0},'o-','MarkerSize',12)
end
hold off
legend('N = 20','N=50','N=100')
legend boxoff
xlabel('input sparsity')
ylabel('sparsity of w')

fileNameEps = [sFolder,filesep,'Law_spW_spc_diffN_gcmi.eps'];
fileNameFig = [sFolder,filesep,'powerLaw_spW_spc_diffN_gcmi.fig'];
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)

% plot the histogram of elments
figure 
hold on
for i0 = 1:length(uniqOdor)
    plot(sparsity{i0},spOut{i0},'o-','MarkerSize',12)
end


% ===================   plot how results depend on alpha ============

str_marker = '2018-06-29';       %folder with this string contains the data we need
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(files,str,'once'));
%for the amplitude data
targetNames = files(FIND(str_marker));


s1 = '(?<= *N)[\d.]+(?=_)';
s2 = '(?<= *_R)[\d]+(?=_)';
s3 = '(?<= *_S)[\d]+(?=_)';
s5 = '(?<= *sig)[\d.]+(?=_)';

numOdor = zeros(length(targetNames),1);
numRecp = zeros(length(targetNames),1);
spInput = zeros(length(targetNames),1);
sig = zeros(length(targetNames),1);
spW = zeros(length(targetNames),1);

figure
for i0 = 1:length(targetNames)
    load(fullfile(dataFolder,targetNames{i0}))
    numOdor(i0) = str2num(char(regexp(targetNames{i0},s1,'match')));
    numRecp(i0) = str2num(char(regexp(targetNames{i0},s2,'match')));
    spInput(i0) = str2num(char(regexp(targetNames{i0},s3,'match')));
    sig(i0) = str2num(char(regexp(targetNames{i0},s5,'match')));
    
    t = allMat(:);
    subplot(3,3,i0)
    
%     histogram(t,'Normalization','probability')
    histogram(t(abs(t)<20),'Normalization','probability')

    lg = ['N=',num2str(numOdor(i0)),', R = ',num2str(numRecp(i0)),', sp = ',num2str(spInput(i0)),...
        ', \sigma = ',num2str(sig(i0))];
    lgh = legend(lg,'Interpreter','latex');
    set(lgh,'FontSize',14)
    legend boxoff
    xlim([-20,5])
    set(gca,'FontSize',16)
    xlabel('$\ln(w)$','Interpreter','latex','FontSize',20)
    ylabel('probability','FontSize',20)
    spW(i0) = sum(t>-15)/length(t);
end

% group the data
uniqAlp = [1.5,2,3];
uniqSp = sort(uniq(spInput));
sparsity = cell(3,1);
spOut = cell(3,1);
for i0 = 1:3
    inx = sig == uniqAlp(i0);
    sparsity{i0} = spInput(inx);
    spOut{i0} = spW(inx);
end

% plot the histogram of elments
figure 
hold on
for i0 = 1:length(uniqAlp)
    plot(sparsity{i0},spOut{i0},'o-','MarkerSize',12)
end
legend('\alpha = 1.5','\alpha = 2','\alpha = 3')
legend boxoff
xlabel('input sparsity')
ylabel('sparsity of w')

%% This part summarize data for the gcmi with different N

% dataFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_Ndp/M20diffN_1021';
dataFolder = '/Users/shan/Dropbox/olfactionProject/data/fig3/GcmiNdpNew1023';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
files = {allFile.name}';

% allN = [30 40 50 80 100 120 150 200];
% allN = [40:10:100,120];
allN = [40 50 60 70 90 100 120 150];  % all the number of odorants
M = 13;
sp = 3;
sig = 2;

s1 = '(?<= *N)[\d.]+(?=_)';
s2 = '(?<= *_R)[\d]+(?=_)';
s3 = '(?<= *_S)[\d]+(?=_)';
s5 = '(?<= *sig)[\d.]+(?=_)';

numOdor = zeros(length(files),1);
numRecp = zeros(length(files),1);
spInput = zeros(length(files),1);
% sig = zeros(length(files),1);
% spW = zeros(length(files),1);
meanFmin = zeros(length(files),1);
stdFmin = zeros(length(files),1);
meanW = zeros(length(files),1);
stdW = zeros(length(files),1);
meanSigW = zeros(length(files),1);
stdSigW = zeros(length(files),1);
meanSpW = zeros(length(files),1);
stdSpW = zeros(length(files),1);

% figure
for i0 = 1:length(files)
    load(fullfile(dataFolder,files{i0}))
    numOdor(i0) = str2num(char(regexp(files{i0},s1,'match')));
    inx = find(allN == numOdor(i0)); % the ordered index

    numRecp(inx) = str2num(char(regexp(files{i0},s2,'match')));
    spInput(inx) = str2num(char(regexp(files{i0},s3,'match')));
%     sig(inx) = str2num(char(regexp(files{i0},s5,'match')));
    
%     meanFmin(inx) = mean(-allfmin);
%     stdFmin(inx) = std(-allfmin);
    
    allActW = nan(size(allMat,2),1);
    allSpW = nan(size(allMat,2),1);
    allSigW = nan(size(allMat,2),1);
    fminSel = nan(size(allMat,2),1);
    
    LB = -4*sig; %lower limit
    UB = 5;  %upper limit
    Fthd = 10; % to determine the simualtion quality
        
    for j0 = 1:size(allMat,2)
%         allActW(j0,1) = mean(temp(temp > -4*sig(inx)));
%         allSigW(j0,1) = std(temp(temp > -4*sig(inx)));
        if allfmin(j0)< Fthd
            fminSel(j0) = allfmin(j0);
            temp = allMat(:,j0);
            t = temp(temp > LB);
            actData = t(t<=UB);
            pd = fitdist(actData,'normal');
            allActW(j0) = pd.mean;
            allSigW(j0) = pd.sigma;
            allSpW(j0) = sum(temp > LB)/length(temp);
        end
    end
    
    meanFmin(inx) = nanmean(fminSel);
    stdFmin(inx) = nanstd(fminSel);
    
    meanW(inx) = nanmean(allActW);
    stdW(inx) = nanstd(allActW);
    
    meanSigW(inx) = nanmean(allSigW);
    stdSigW(inx) = nanstd(allSigW);
    
    meanSpW(inx) = nanmean(allSpW);
    stdSpW(inx) = nanstd(allSpW);
%     meanSpW(inx) = nanmean(sum(allMat > thd,1))/size(allMat,1);
%     stdSpW(inx) = nanstd(sum(allMat > thd,1)/size(allMat,1));
    
end

% save the data
% put all variables into a struct
dataSumm = struct('allN',allN,'meanFmin',meanFmin,'stdFmin',stdFmin,'meanW',meanW,...
    'stdW',stdW,'meanSigW',meanSigW,'stdSigW',stdSigW,'meanSpW',meanSpW,'stdSpW',...
    stdSpW,'nRecp',numRecp,'sp',spInput,'sig',sig);
dName = fullfile(sFolder,['gcmi_Ndp_','M',num2str(M),'Sp',num2str(sp),...
    'sig',num2str(sig),'_',date,'.mat']);
save(dName,'dataSumm')

%% This part summarize data for the gcmi with different M

dataFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_Mdp';
% dataFolder = '/Users/shan/Dropbox/olfactionProject/data/GcmiMdp12-09';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*2018-12-24.mat'));
files = {allFile.name}';

N = 50;
% allM = [10,20:5:50];  % all the number of odorants
allM = 10:5:40;
sp = 5;    % this need to be modified accordingly
sig = 2;

s1 = '(?<= *N)[\d.]+(?=_)';
s2 = '(?<= *_R)[\d]+(?=_)';
s3 = '(?<= *_S)[\d]+(?=_)';
s5 = '(?<= *sig)[\d.]+(?=_)';

numOr = zeros(length(files),1);
numRecp = zeros(length(files),1);
spInput = zeros(length(files),1);

meanFmin = zeros(length(files),1);
stdFmin = zeros(length(files),1);
meanW = zeros(length(files),1);
stdW = zeros(length(files),1);
meanSigW = zeros(length(files),1);
stdSigW = zeros(length(files),1);
meanSpW = zeros(length(files),1);
stdSpW = zeros(length(files),1);

lb = [-4*sig*ones(1,5),-5*sig,-6*sig];
ub = [2,2,2,5,5,5,6];

% figure
for i0 = 1:length(files)
    load(fullfile(dataFolder,files{i0}))
    numOr(i0) = str2num(char(regexp(files{i0},s2,'match')));
    inx = find(allM == numOr(i0)); % the ordered index

    numRecp(inx) = str2num(char(regexp(files{i0},s2,'match')));
    spInput(inx) = str2num(char(regexp(files{i0},s3,'match')));
%     sig(inx) = str2num(char(regexp(files{i0},s5,'match')));
    
%     meanFmin(inx) = mean(-allfmin);
%     stdFmin(inx) = std(-allfmin);
    
    allActW = nan(size(allMat,2),1);
    allSpW = nan(size(allMat,2),1);
    allSigW = nan(size(allMat,2),1);
    fminSel = nan(size(allMat,2),1);
    
%     LB = -4*sig;    %lower limit
%     UB = 5;         %upper limit
    LB = lb(i0);    %lower limit
    UB = ub(i0);         %upper limit
    Fthd = 50;      % to determine the simualtion quality
        
    for j0 = 1:size(allMat,2)
%         allActW(j0,1) = mean(temp(temp > -4*sig(inx)));
%         allSigW(j0,1) = std(temp(temp > -4*sig(inx)));
        if allfmin(j0)< Fthd
            fminSel(j0) = allfmin(j0);
            temp = allMat(:,j0);
            t = temp(temp > LB);
            actData = t(t<=UB);
            pd = fitdist(actData,'normal');
            allActW(j0) = pd.mean;
            allSigW(j0) = pd.sigma;
            allSpW(j0) = sum(temp > LB)/length(temp);
        end
    end
    
    meanFmin(inx) = nanmean(fminSel);
    stdFmin(inx) = nanstd(fminSel);
    
    meanW(inx) = nanmean(allActW);
    stdW(inx) = nanstd(allActW);
    
    meanSigW(inx) = nanmean(allSigW);
    stdSigW(inx) = nanstd(allSigW);
    
    meanSpW(inx) = nanmean(allSpW);
    stdSpW(inx) = nanstd(allSpW);
%     meanSpW(inx) = nanmean(sum(allMat > thd,1))/size(allMat,1);
%     stdSpW(inx) = nanstd(sum(allMat > thd,1)/size(allMat,1));
    
end

% save the data
% put all variables into a struct
dataSumm = struct('allM',allM,'meanFmin',meanFmin,'stdFmin',stdFmin,'meanW',meanW,...
    'stdW',stdW,'meanSigW',meanSigW,'stdSigW',stdSigW,'meanSpW',meanSpW,'stdSpW',...
    stdSpW,'nRecp',numRecp,'sp',spInput,'sig',sig);
dName = fullfile(sFolder,['gcmi_Mdp_','N',num2str(N),'Sp',num2str(sp),...
    'sig',num2str(sig),'_',date,'.mat']);
save(dName,'dataSumm')


%% This part summarize data for gcmi with differnt sigma_c
% dataFolder = '/Users/shan/Dropbox/olfactionProject/data/fig3/GcmiDifcons10-15';
dataFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_sig';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
filesRaw = {allFile.name}';

% screen the files, only part are used. Since the folder may contain
% heterogeneous data set
str_marker = '2018-10-21';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

% define the basic parameters
% R = [60,100,150,200];
num = 40;    % number of repeat
M = 20;      % number of receptors
N = 100;     % number of odorants
sp = 3;
% allSig = 1.1:0.3:5.9;
allSig = [1.1, 1.5:0.5:5];


s1 = '(?<= *N)[\d.]+(?=_)';
s2 = '(?<= *_R)[\d]+(?=_)';
s3 = '(?<= *_S)[\d]+(?=_)';
s5 = '(?<= *sig)[\d.]+(?=_)';

% allSig = zeros(length(files),1);
numRecp = zeros(length(files),1);
spInput = zeros(length(files),1);
% sig = zeros(length(files),1);
% spW = zeros(length(files),1);
meanFmin = zeros(length(files),1);
stdFmin = zeros(length(files),1);
meanW = zeros(length(files),1);
stdW = zeros(length(files),1);
meanSigW = zeros(length(files),1);
stdSigW = zeros(length(files),1);
meanSpW = zeros(length(files),1);
stdSpW = zeros(length(files),1);

% figure
for i0 = 1:length(files)
    load(fullfile(dataFolder,files{i0}))
    sig = str2num(char(regexp(files{i0},s5,'match')));
    inx = find(round(10*allSig) == round(10*sig)); % the ordered index

%     numRecp(inx) = str2num(char(regexp(files{i0},s2,'match')));
%     spInput(inx) = str2num(char(regexp(files{i0},s3,'match')));
%     sig(inx) = str2num(char(regexp(files{i0},s5,'match')));
    
    
    allActW = nan(size(allMat,2),1);
    allSpW = nan(size(allMat,2),1);
    allSigW = nan(size(allMat,2),1);
    fminSel = nan(size(allMat,2),1);
    
    LB = -4*sig; %lower limit
    UB = 5;  %upper limit
    Fthd = 30; % to determine the simualtion quality, this differs from situations
        
    for j0 = 1:size(allMat,2)

        if allfmin(j0)< Fthd
            fminSel(j0) = allfmin(j0);
            temp = allMat(:,j0);
            t = temp(temp > LB);
            actData = t(t<=UB);
            pd = fitdist(actData,'normal');
            allActW(j0) = pd.mean;
            allSigW(j0) = pd.sigma;
            allSpW(j0) = sum(temp > LB)/length(temp);
        end
    end
    
    meanFmin(inx) = nanmean(fminSel);
    stdFmin(inx) = nanstd(fminSel);
    
    meanW(inx) = nanmean(allActW);
    stdW(inx) = nanstd(allActW);
    
    meanSigW(inx) = nanmean(allSigW);
    stdSigW(inx) = nanstd(allSigW);
    
    meanSpW(inx) = nanmean(allSpW);
    stdSpW(inx) = nanstd(allSpW);
    
end

% save the data
% put all variables into a struct
dataSumm = struct('meanFmin',meanFmin,'stdFmin',stdFmin,'meanW',meanW,...
    'stdW',stdW,'meanSigW',meanSigW,'stdSigW',stdSigW,'meanSpW',meanSpW,'stdSpW',...
    stdSpW,'nRecp',M,'sp',sp,'allSig',allSig);

dName = fullfile(saveFolder,['gcmi_sigdp_','N',num2str(N),'M',num2str(M),'Sp',num2str(sp),...
    '_',date,'.mat']);
save(dName,'dataSumm')


% data structure store all the data


%% This part summarize data for gcmi with differnt input sparsity n
dataFolder = '/Users/shan/Dropbox/olfactionProject/data/fig3/GcmiSdp10-17';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
filesRaw = {allFile.name}';

% screen the files, only part are used. Since the folder may contain
% heterogeneous data set
str_marker = '2018-10-17';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

% define the basic parameters
% R = [60,100,150,200];
num = 40;    % number of repeat
M = 13;      % number of receptors
N = 50;     % number of odorants
sig = 2;
allSp = 2:1:10; % all different input sparsity

s1 = '(?<= *N)[\d.]+(?=_)';
s2 = '(?<= *_R)[\d]+(?=_)';
s3 = '(?<= *_S)[\d]+(?=_)';
s5 = '(?<= *sig)[\d.]+(?=_)';

% allSig = zeros(length(files),1);
% numRecp = zeros(length(files),1);
% spInput = zeros(length(files),1);
% sig = zeros(length(files),1);
% spW = zeros(length(files),1);
meanFmin = zeros(length(files),1);
stdFmin = zeros(length(files),1);
meanW = zeros(length(files),1);
stdW = zeros(length(files),1);
meanSigW = zeros(length(files),1);
stdSigW = zeros(length(files),1);
meanSpW = zeros(length(files),1);
stdSpW = zeros(length(files),1);

% figure
for i0 = 1:length(files)
    load(fullfile(dataFolder,files{i0}))
    sp = str2num(char(regexp(files{i0},s3,'match')));
    inx = find(allSp == sp); % the ordered index

%     numRecp(inx) = str2num(char(regexp(files{i0},s2,'match')));
%     spInput(inx) = str2num(char(regexp(files{i0},s3,'match')));
%     sig(inx) = str2num(char(regexp(files{i0},s5,'match')));
    
    
    allActW = nan(size(allMat,2),1);
    allSpW = nan(size(allMat,2),1);
    allSigW = nan(size(allMat,2),1);
    fminSel = nan(size(allMat,2),1);
    
    LB = -4*sig; %lower limit
    UB = 10;  %upper limit
    Fthd = 4; % to determine the simualtion quality, this differs from situations
        
    for j0 = 1:size(allMat,2)

        if allfmin(j0)< Fthd
            fminSel(j0) = allfmin(j0);
            temp = allMat(:,j0);
            t = temp(temp > LB);
            actData = t(t<=UB);
            pd = fitdist(actData,'normal');
            allActW(j0) = pd.mean;
            allSigW(j0) = pd.sigma;
            allSpW(j0) = sum(temp > LB)/length(temp);
        end
    end
    
    meanFmin(inx) = nanmean(fminSel);
    stdFmin(inx) = nanstd(fminSel);
    
    meanW(inx) = nanmean(allActW);
    stdW(inx) = nanstd(allActW);
    
    meanSigW(inx) = nanmean(allSigW);
    stdSigW(inx) = nanstd(allSigW);
    
    meanSpW(inx) = nanmean(allSpW);
    stdSpW(inx) = nanstd(allSpW);
    
end

% save the data
% put all variables into a struct
dataSumm = struct('meanFmin',meanFmin,'stdFmin',stdFmin,'meanW',meanW,...
    'stdW',stdW,'meanSigW',meanSigW,'stdSigW',stdSigW,'meanSpW',meanSpW,'stdSpW',...
    stdSpW,'nRecp',M,'sp',allSp,'sig',sig);

dName = fullfile(saveFolder,['gcmi_spWdp_','N',num2str(N),'M',num2str(M),'Sp',num2str(sp),...
    '_',date,'.mat']);
save(dName,'dataSumm')


% data structure store all the data


%% This part summarize gcmi whole matrix simulation data, from Qianyi
dataFolder = '/Users/shan/Dropbox/olfactionProject/data/GcmiDifSigData';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
filesRaw = {allFile.name}';

% screen the files, only part are used. Since the folder may contain
% heterogeneous data set
str_marker = '2018-05-18';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

% define the basic parameters
allN = [10 20 30 50 100 150 200];
R = 30;      % number of receptors
num = 20;    % number of repeat
sp = 2;      % sparsity of input
allSig = 1.1:0.3:3.2;  % different sigma

for i0 = 1:length(files)
    s1 = '(?<= *N)[\d.]+(?=_)';
    s3 = '(?<= *_sig)[\d.]+(?=_)';
    N = str2num(char(regexp(files{i0},s1,'match')));
    sig = str2num(char(regexp(files{i0},s3,'match')));
    inxN = find(allN == N);
    inxSig = find(round(allSig*10) == round(sig*10));
    
    load(char(fullfile(dataFolder,filesep,files{i0})));
    if num ~= length(allfmin)
        disp('number of repeats in this simulation do not match!')
    end
    
    %hard threshold of sparisty
    thd = -5*sig;
    inx = allfmin < 2e2;  % only preserve the successed!
    allMeanFmin(inxSig,inxN) = mean(-allfmin(inx));
    allStdFmin(inxSig,inxN) = std(-allfmin(inx));
      
    meanW = [];
    stdW = [];
    for k0 = 1 :length(inx)
        temp = allMat(:,inx(k0));
        meanW = [meanW,mean(temp(temp > thd))];
        stdW = [stdW,std(temp(temp > thd))];     
    end
    allMeanAveW(inxSig,inxN) = mean(meanW);
    allMeanSigW(inxSig,inxN) = mean(stdW);
    allStdAveW(inxSig,inxN) = std(meanW);
    allStdSigW(inxSig,inxN) = std(stdW);
    
    temp2 = allMat(:,inx);
    allMeanSpW(inxSig,inxN) = mean(sum(temp2>thd,1)/size(allMat,1));
    allStdSpW(inxSig,inxN) = std(sum(temp2>thd,1)/size(allMat,1));

end
% save the data
saveName = ['gcmi_diffN_diffSig_M',num2str(R),'Summary_QY_',date,'.mat'];
dataName = fullfile(saveFolder,saveName);
save(dataName,'R','allN','sp','allSig','allMeanFmin','allStdFmin','allMeanAveW',...
    'allStdAveW','allStdSigW','allMeanSpW','allStdSpW')

%% This part analyzes the data with both inhibition and excitation
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_inhi/gcmiInhiN50M10diffbasal_1007';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
allFile = dir([dFolder,filesep,'*.mat']);
allfiles = {allFile.name}';

str_marker = '2018-10-07';   
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(allfiles,str,'once'));
files = allfiles(FIND(str_marker));

N = 50;
M = 10;
sp = 3;
% sig = 2;
r0 = 0.05;

% default graphics settings
defaultGraphicsSetttings
allRatio = zeros(length(files),1);
stdRatio = zeros(length(files),1);
meanFmin = zeros(length(files),1);
stdFmin = zeros(length(files),1);
% allr0 = zeros(length(files),1);
exciEle = zeros(length(files),2); %mean and std
inhiEle = zeros(length(files),2); %mean and std
allSig = zeros(length(files),1);

% set the color
lpurple = [230,61,144]/256;
lblue = [46,167,224]/256;
myGreen = [121,191,132]/256;
myRed = [234,85,65]/256;
% figure
figFlag = 0;  % plot or not
for i0 = 1:length(files)
    str = '(?<= *sig)[\d.]+(?=_)';
    allSig(i0) = str2num(char(regexp(files{i0},str,'match')));
    load(files{i0})
    
    allRatio(i0) = sum(allSign(:)  < 0)/length(allMat(:));
    stdRatio(i0) = std(sum(allSign < 0,1)/length(allMat(:,1)));
    meanFmin(i0) = - mean(allfmin(allfmin < 2));
    stdFmin(i0) = std(allfmin(allfmin < 10));
    exciEle(i0,:) =  [mean(log(allMat(allSign > 0))),std(log(allMat(allSign > 0)))];
    inhiEle(i0,:) =  [mean(log(allMat(allSign < 0))),std(log(allMat(allSign < 0)))];
%     histogram(log(allMat(allSign > 0)),'Normalization','probability')
%     histogram(log(allMat(allSign  < 0)),'Normalization','probability')
%     subplot(3,4,i0)

if figFlag
    figure
    histogram(log(allMat(allSign > 0)),'EdgeColor','none','FaceColor',lpurple)
    hold on
    histogram(log(allMat(allSign  < 0)),'EdgeColor','none','FaceColor',lblue)
    
    % add a patch
    xp1 = [3.8,4.5,4.5,3.8];
    yp1 = [600,600,900,900];
    v1 = [3.5, 600;4.2, 600;4.2, 900;3.5, 900];
    f1 = [1 2 3 4];
    patch('Faces',f1,'Vertices',v1,...
    'EdgeColor','k','FaceColor','none','LineWidth',2);
    
    v2 = [3.52, 600;4.18, 600;4.18, 600 + 300*(allr0(i0));3.52, 600 + 300*(allr0(i0))];
    f2 = [1 2 3 4];
    patch('Faces',f2,'Vertices',v2,...
    'EdgeColor','none','FaceColor',myRed);
%     patch(xp1,yp1,myGreen)
    hold off
    xp = 2.2;
    yp = 940;
    txt = ['$r_0$ =', num2str(allr0(i0))];
    text(xp,yp,txt,'Interpreter','latex','FontSize',20)

    xlim([-6 5])
    ylim([0,1000])
%     if i0==1
%         legend('excitation','inhibition')
%         legend boxoff
%     end
    lg = legend('excitation','inhibition','Location','northwest');
    set(lg,'FontSize',16)
    legend boxoff
    xlabel('$\ln(w)$','Interpreter','latex')
    ylabel('counts')
    prefix = ['inhi_acti_histW_',num2str(i0),'.png'];
    print(gcf,'-r600','-dpng',fullfile(sFolder,prefix))
end
end

% plot how fmin changes with spontaneous activity

%summarize and save the data
% [uniqr0,inx] = sort(allr0);  %sort based on the 
[uniqSig,inx] = sort(allSig);  %sort based on the 

summData = struct('N',N,'M',M,'sp',sp,'allSig',uniqSig,'meanFmin',...
    meanFmin(inx),'stdFmin',stdFmin(inx),'allRatio',allRatio(inx),'stdRatio',...
    stdRatio(inx),'exciEle',exciEle(inx,:),'inhiEle',inhiEle(inx,:));
dFile = fullfile(sFolder, ['gcmiInhi_N',num2str(N),'_M',num2str(M),'_sp',...
    num2str(sp),'_r0',num2str(r0),'_',date,'.mat']);
save(dFile,'summData')

% allr0 = [0.1,0.15,0.2:0.1:0.8,0.85,0.9];
figure
errorbar(uniqSig',meanFmin(inx),stdFmin(inx),'o-','MarkerSize',12,'LineWidth',2)
xlabel('$\sigma_c$','Interpreter','latex')
ylabel('$f_{min}$','Interpreter','latex')
ylim([-2,-1])

fileNameEps = fullfile(sFolder,['gcmiInhi_N',num2str(N),'_M',num2str(M),'_sp',...
    num2str(sp),'_sig',num2str(sig),'_inhi_fmin_',date,'.eps']);
fileNameFig = fullfile(sFolder,['gcmiInhi_N',num2str(N),'_M',num2str(M),'_sp',...
    num2str(sp),'_sig',num2str(sig),'_inhi_fmin_',date,'.fig']);
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)


% how inhibitory ratio changes with r_0
figure
hold on
errorbar(uniqSig,allRatio(inx),stdRatio(inx),'o-','MarkerSize',12,'LineWidth',2)
% plot([0,1],[0,1],'k--','LineWidth',1)
hold off
ylim([0.1 0.25])
xlabel('$\sigma_c$','Interpreter','latex')
ylabel('$\rho$','Interpreter','latex')
fileNameEps = fullfile(sFolder,['gcmiInhi_N',num2str(N),'_M',num2str(M),'_sp',...
    num2str(sp),'_sig',num2str(sig),'_inhi_ratio_',date,'.eps']);
fileNameFig = fullfile(sFolder,['gcmiInhi_N',num2str(N),'_M',num2str(M),'_sp',...
    num2str(sp),'_sig',num2str(sig),'_inhi_ratio_',date,'.fig']);
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)


% distribution of excitatory and inhibitory elements
% firs the mean value
figure
hold on
plot(uniqSig,exciEle(inx,1),'o-','MarkerSize',12,'LineWidth',2)
plot(uniqSig,inhiEle(inx,1),'o-','MarkerSize',12,'LineWidth',2)
hold off
legend('exci','inhi')
legend boxoff
xlabel('$\sigma_c$','Interpreter','latex')
ylabel('$\mu_w$','Interpreter','latex')
fileNameEps = fullfile(sFolder,['gcmiInhi_N',num2str(N),'_M',num2str(M),'_sp',...
    num2str(sp),'_sig',num2str(sig),'_inhi_meanVal_',date,'.eps']);
fileNameFig = fullfile(sFolder,['gcmiInhi_N',num2str(N),'_M',num2str(M),'_sp',...
    num2str(sp),'_sig',num2str(sig),'_inhi_meanVal_',date,'.fig']);
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)

figure
hold on
plot(uniqSig,exciEle(inx,2),'o-','MarkerSize',12,'LineWidth',2)
plot(uniqSig,inhiEle(inx,2),'o-','MarkerSize',12,'LineWidth',2)
hold off
legend('exci','inhi')
legend boxoff
ylim([0,2])
xlabel('$\sigma_c$','Interpreter','latex')
ylabel('$\sigma_w$','Interpreter','latex')

fileNameEps = fullfile(sFolder,['gcmiInhi_N',num2str(N),'_M',num2str(M),'_sp',...
    num2str(sp),'_sig',num2str(sig),'_inhi_sigW_',date,'.eps']);
fileNameFig = fullfile(sFolder,['gcmiInhi_N',num2str(N),'_M',num2str(M),'_sp',...
    num2str(sp),'_sig',num2str(sig),'_inhi_sigW_',date,'.fig']);
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)


% video
file1 = 'int_N100_R10_S5_sig2_alp0.1_2018-06-14.mat';
file2 = 'int_N100_R10_S5_sig2_alp0.2_2018-06-14.mat';
load(fullfile(dFolder,file2))
rho1 =  sum(allSign(:)  < 0)/length(allMat(:));
fmin1 = - mean(allfmin(allfmin < 10));
figure
hold on
histogram(log(allMat(allSign > 0)))
histogram(log(allMat(allSign < 0)))
hold off
title(['f_{min}=',num2str(fmin1),',\rho = ',num2str(rho1)],'FontSize',20)
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('counts')

alp = [1.5,2,3];
c0 = 1;
X = 10.^(0:0.02:5);

Y1 = c0*(X/c0).^(1-alp(1));
Y2 = c0*(X/c0).^(1-alp(2));
Y3 = c0*(X/c0).^(1-alp(3));

figure
hold on
plot(X,Y1)
plot(X,Y2)
plot(X,Y3)
hold off
legend('\alpha = 1.5','\alpha = 2','\alpha = 3')
legend boxoff
set(gca,'XScale','log','YScale','log')
xlabel('$\ln(c)$','Interpreter','latex')
ylabel('$\ln(p)$','Interpreter','latex')


%% This part analyzes the data using gcmi with distribution 
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_distr';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
allFile = dir([dFolder,filesep,'*.mat']);
allfiles = {allFile.name}';

str_marker = '2018-08-28';   
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(allfiles,str,'once'));
files = allfiles(FIND(str_marker));

% default graphics settings
allN = 100;
allSp = 2;
% allSig = [1.1,1.5,2.1,2.7,3.1,3.5];%08-23-2018
allSig = [1.1,1.5,2.1,2.5,2.9,3.3,3.7,4.1];%08-23-2018


meanFmin = zeros(length(files),1);
stdFmin = zeros(length(files),1);
meanRho = zeros(length(files),1);
stdRho = zeros(length(files),1);
meanW = zeros(length(files),1);
stdW = zeros(length(files),1);
meanSigW= zeros(length(files),1);
stdSigW = zeros(length(files),1);

% meanFmin = zeros(length(allN),length(allSp));
% stdFmin = zeros(length(allN),length(allSp));
% meanRho = zeros(length(allN),length(allSp));
% stdRho = zeros(length(allN),length(allSp));
% meanW = zeros(length(allN),length(allSp));
% stdW = zeros(length(allN),length(allSp));
% meanSigW = zeros(length(allN),length(allSp));
% stdSigW = zeros(length(allN),length(allSp));

% allSp = zeros(length(files),1);

for i0 = 1:length(files)
    str = '(?<= *_S)[\d.]+(?=_)';
%     str2 = '(?<= *_N)[\d.]+(?=_)';
    str3 = '(?<= *_sig)[\d.]+(?=_)';
    sig = str2num(char(regexp(files{i0},str3,'match')));
%     nOdor = str2num(char(regexp(files{i0},str2,'match')));
%     inx1 = find(allN == nOdor);
    inx2 = 1;
    inx1 = find(allSig == sig);
    load(files{i0})
    inx = allfmin < 1e3;  % only preserve the successed!
    meanFmin(inx1,inx2) = mean(-allfmin(inx));
    stdFmin(inx1,inx2) = std(-allfmin(inx));
    temp = mean(allParam(:,inx),2);
    meanW(inx1,inx2) = temp(1);
    meanSigW(inx1,inx2) = temp(2);
    meanRho(inx1,inx2) = temp(3);
    temp2 = std(allParam(:,inx),0,2);
    stdW(inx1,inx2) = temp2(1);
    stdSigW(inx1,inx2) = temp2(2);
    stdRho(inx1,inx2) = temp2(3);
end

% reorder the data
summData=struct('allN',allN','allSp',allSp','allSig',allSig,'meanFmin',meanFmin,...
    'stdFmin',stdFmin,'meanW',meanW,'stdW',stdW,'meanSigW',meanSigW,'stdSigW',stdSigW,...
    'meanRho',meanRho,'stdRho',stdRho);
dName=['gcmi_distri_summData_',date,'.mat'];
% save(fullfile(sFolder,dName),'-struct','summData');
save(fullfile(sFolder,dName),'summData');

% uniqN = sort(uniq(allN));
% figure
% inx = allSp ==2;
% plot(N',-meanFmin(inx))


%% this part analyze data for partial inhibitory receptors
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_partInhi';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
allFile = dir([dFolder,filesep,'*.mat']);
allfiles = {allFile.name}';

str_marker = '2019-06-25';   
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(allfiles,str,'once'));
files = allfiles(FIND(str_marker));

% default graphics settings
allN = 50;
M = 20;
allSp = 3;
sig = 2;
% r0 = [0.2,0.5];     % basal activity
r0 = 0.4;     % maximum basal activity
frac = 0.75;   % fraction of inhibtory receptors


allWi = cell(length(r0),length(frac));
Sign = cell(length(r0),length(frac));
allWe = cell(length(r0),length(frac));

meanFmin = zeros(length(r0),length(frac));
stdFmin = zeros(length(r0),length(frac));
meanRho = zeros(length(r0),length(frac));
stdRho = zeros(length(r0),length(frac));
meanW = zeros(length(r0),length(frac));
stdW = zeros(length(r0),length(frac));
meanWe = zeros(length(r0),length(frac));
stdWe = zeros(length(r0),length(frac));
meanWi = zeros(length(r0),length(frac));
stdWi = zeros(length(r0),length(frac));
meanWie = zeros(length(r0),length(frac)); %excitatory in first part
stdWie = zeros(length(r0),length(frac));

meanSigW = zeros(length(r0),length(frac));
stdSigW = zeros(length(r0),length(frac));

meanSigWe = zeros(length(r0),length(frac));
stdSigWe = zeros(length(r0),length(frac));

meanSigWi = zeros(length(r0),length(frac));
stdSigWi = zeros(length(r0),length(frac));
meanSigWie = zeros(length(r0),length(frac));  % standard devitation of excitatory in first part
stdSigWie = zeros(length(r0),length(frac));

% allSp = zeros(length(files),1);

for i0 = 1:length(files)
    str1 = '(?<= *_alp)[\d.]+(?=_)';
    str2 = '(?<= *_frac)[\d.]+(?=_)';
    alp = str2num(char(regexp(files{i0},str1,'match')));
    fracInhi = str2num(char(regexp(files{i0},str2,'match')));
    inx1 = find(r0 == alp);
    inx2 = find(frac == fracInhi);
    load(fullfile(dFolder,files{i0}))
    
    inx = allfmin < 10;  % only preserve the successed!
    meanFmin(inx1,inx2) = mean(-allfmin(inx));
    stdFmin(inx1,inx2) = std(-allfmin(inx));
    
    allWi{inx1,inx2} = allInhiW;
    Sign{inx1,inx2} = allSign;
    allWe{inx1,inx2} = allExciW;
    
    % mean and std of excitatory part
    meanWe(inx1,inx2) = mean(log(allExciW(:)));
    stdWe(inx1,inx2) = std(mean(log(allExciW),1));
    meanWi(inx1,inx2) = mean(log(allInhiW(allSign<0)));
    temp = nan(size(allInhiW,1),size(allInhiW,2));
    temp(allSign < 0) = log(allInhiW(allSign<0));
    stdWi(inx1,inx2) = nanstd(nanmean(temp,1));
    
    meanSigWe(inx1,inx2) = mean(std(log(allExciW),0,1));
    stdSigWe(inx1,inx2) = std(std(log(allExciW),0,1));
    meanSigWi(inx1,inx2) = mean(nanstd(temp,0,1));
    stdSigWi(inx1,inx2) = std(nanstd(temp,0,1));
    
    temp2 = nan(size(allInhiW,1),size(allInhiW,1));
    temp2(allSign > 0) = log(allInhiW(allSign>0));
    meanWie(inx1,inx2) = mean(log(allInhiW(allSign>0)));
    stdWie(inx1,inx2) = nanstd(nanmean(temp2,1));
    meanSigWie(inx1,inx2) = nanmean(nanstd(temp2,0,1));
    stdSigWie(inx1,inx2) = std(nanstd(temp2,1));
    

    meanRho(inx1,inx2) = mean(sum(allSign < 0,1)/size(allSign,1)); 
    stdRho(inx1,inx2) = std(sum(allSign < 0,1)/size(allSign,1));
    

end

% reorder the data
summData=struct('allr0',r0','allfrac',frac','meanFmin',meanFmin,...
    'stdFmin',stdFmin,'meanWe',meanWe,'stdWe',stdWe,'meanWi',meanWi,'stdWi',...
    stdWi,'meanWie',meanWie,'stdWie',stdWie,'meanSigWe',meanSigWe,'stdSigWe',...
    stdSigWe,'meanSigWi',meanSigWi,'meanSigWie',meanSigWie,'meanRho',meanRho);
%'allWi',allWi,'allWe',allWe,'allWie',allWie
dName=['gcmi_partInhi_summData_N',num2str(allN),'M',num2str(M),'sp',num2str(allSp),'sig',num2str(sig),'_',date,'.mat'];
save(fullfile(sFolder,dName),'allWi','allWe','Sign','summData');
% uniqN = sort(uniq(allN));
% figure
% inx = allSp ==2;
% plot(N',-meanFmin(inx))

%% This part summize the data for partial inhibition with different spontaneous firing rate
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_partInhi';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

file = 'gcmiPartInhi_N50_R20_S2_sig2_alp0.4_frac0.75_2019-06-25.mat';

% default graphics settings
allN = 50;
M = 20;
allSp = 2;
sig = 2;
r0 = 0.4;     % maximum basal activity
frac = 0.75;   % fraction of inhibtory receptors

load(fullfile(dFolder,file))

% check the qulity of the data
if sum(allfmin > 10)>0
   disp('Causion! the quality of data is low!')
end
meanFmin = mean(allfmin);
stdFmin = std(allfmin); 

repeat = size(allSign,2);  % number of repeats
temp = reshape(allSign,ceil(M*frac),allN,repeat);
meanInhiRatio = mean(sum(temp<0,2),3)/allN;  % average number of inhibiton ratio
stdInhiRatio = std(sum(temp<0,2)/allN,0,3);  % standard deviation of inhibitory ratio

% calculate the non-zero elements of the excitatory part matrix
thd = -10;  % a hard threshold to specify the zero elements in the excitatory part matri
menaRho = sum(log(allExciW(:)) < thd)/length(allExciW(:)); % mean zero element fraction
stdRho = std(sum(log(allExciW)<-10)/size(allExciW,1));  % standard deviation


% reorder the data
summData=struct('meanInhiRatio',meanInhiRatio,'stdInhiRatio',stdInhiRatio,...
    'menaRho',menaRho,'stdRho',stdRho);
%'allWi',allWi,'allWe',allWe,'allWie',allWie
dName=['gcmi_partInhi_summData_N',num2str(allN),'M',num2str(M),'sp',num2str(allSp),...
    'sig',num2str(sig),'frac_',num2str(frac),'_',date,'.mat'];
save(fullfile(sFolder,dName),'allSign','allInhiW','allExciW','summData');


%% this part analyze the data for gcmi distribution
dFolder = '/Users/shan/Dropbox/olfactionProject/data/gcmi/gcmi_distr/diffSig0806';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
allFile = dir([dFolder,filesep,'*.mat']);
allfiles = {allFile.name}';

str_marker = '2018-08-06';   
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(allfiles,str,'once'));
files = allfiles(FIND(str_marker));

% default graphics settings
allN = 100;
% allN = [20 30 40 50 80 100 120 150 200];
allM = 30;
allSp = 5;
% allSig = 2;
allSig = 1.1:0.2:4.1;


% allFmin = cell(length(allM),length(allSp));
% allMeanSpW = cell(length(allM),length(allSp));
% allStdSpW = cell(length(allM),length(allSp));
% 
% allMeanAveW = cell(length(allM),length(allSp));
% allStdAveW = cell(length(allM),length(allSp));

%{
allMeanFmin = cell(length(allN),1);
allStdFmin = cell(length(allN),1);

allMeanSpW = cell(length(allN),1);
allStdSpW = cell(length(allN),1);

allMeanAveW = cell(length(allN),1);
allStdAveW = cell(length(allN),1);

allMeanSigW = cell(length(allN),1);
allStdSigW = cell(length(allN),1);


% initialize
for i0 = 1:length(allN)
    allMeanFmin{i0} = nan(length(allM),length(allSp));
    allStdFmin{i0} = nan(length(allM),length(allSp));
    allMeanSpW{i0} = nan(length(allM),length(allSp));
    allStdSpW{i0} = nan(length(allM),length(allSp));
    allMeanAveW{i0} = nan(length(allM),length(allSp));
    allStdAveW{i0} = nan(length(allM),length(allSp));
    allMeanSigW{i0} = nan(length(allM),length(allSp));
    allStdSigW{i0} = nan(length(allM),length(allSp));
end

% loop through differnt data files
for i0 = 1:length(files)
       
    s1 = '(?<= *N)[\d.]+(?=_)';
    s2 = '(?<= *_R)[\d]+(?=_)';
    s3 = '(?<= *_S)[\d]+(?=_)';
    N = str2num(char(regexp(files{i0},s1,'match')));
    M = str2num(char(regexp(files{i0},s2,'match')));
    sp = str2num(char(regexp(files{i0},s3,'match')));
    ixN = find(allN == N);
    ixM = find(allM == M);
    ixSp = find(allSp == sp);
    
    load(files{i0})
    inx = allfmin < 1e2;  % only preserve the successed!
    allMeanFmin{ixN}(ixM,ixSp) = mean(-allfmin(inx));
    allStdFmin{ixN}(ixM,ixSp) = std(-allfmin(inx));
    
    temp = mean(allParam(:,inx),2);
    allMeanAveW{ixN}(ixM,ixSp) = temp(1);
    allMeanSigW{ixN}(ixM,ixSp) = temp(2);
    allMeanSpW{ixN}(ixM,ixSp) = temp(3);
    
    temp2 = std(allParam(:,inx),0,2);
    allStdAveW{ixN}(ixM,ixSp) = temp2(1);
    allStdSigW{ixN}(ixM,ixSp) = temp2(2);
    allStdSpW{ixN}(ixM,ixSp) = temp2(3);
end

%}



% ======================================================
% change with sigma
% ======================================================
allMeanFmin = zeros(length(allSig),1);
allStdFmin = zeros(length(allSig),1);

allMeanSpW = zeros(length(allSig),1);
allStdSpW = zeros(length(allSig),1);
allMeanAveW = zeros(length(allSig),1);
allStdAveW = zeros(length(allSig),1);

allMeanSigW = zeros(length(allSig),1);
allStdSigW = zeros(length(allSig),1);

for i0 = 1:length(files)
       
    s4 = '(?<= *_sig)[\d.]+(?=_)';
    sig = str2num(char(regexp(files{i0},s4,'match')))
    ixSig = find(round(allSig*10) == round(sig*10))
    
    
    load(files{i0})
    inx = allfmin < 1e2;  % only preserve the successed!
    allMeanFmin(ixSig) = mean(-allfmin(inx));
    allStdFmin(ixSig) = std(-allfmin(inx));
    
    temp = mean(allParam(:,inx),2);
    allMeanAveW(ixSig) = temp(1);
    allMeanSigW(ixSig) = temp(2);
    allMeanSpW(ixSig) = temp(3);
    
    temp2 = std(allParam(:,inx),0,2);
    allStdAveW(ixSig) = temp2(1);
    allStdSigW(ixSig) = temp2(2);
    allStdSpW(ixSig) = temp2(3);
end


%{
% ======================================================
% change with N
% ======================================================

allMeanFmin = zeros(length(allN),1);
allStdFmin = zeros(length(allN),1);

allMeanSpW = zeros(length(allN),1);
allStdSpW = zeros(length(allN),1);
allMeanAveW = zeros(length(allN),1);
allStdAveW = zeros(length(allN),1);

allMeanSigW = zeros(length(allN),1);
allStdSigW = zeros(length(allN),1);

for i0 = 1:length(files)
       
    s1 = '(?<= *_N)[\d]+(?=_)';
    N = str2num(char(regexp(files{i0},s1,'match')))
    ixN = find(allN == N)
    
    
    load(files{i0})
    inx = allfmin < 1e2;  % only preserve the successed!
    allMeanFmin(ixN) = mean(-allfmin(inx));
    allStdFmin(ixN) = std(-allfmin(inx));
    
    temp = mean(allParam(:,inx),2);
    allMeanAveW(ixN) = temp(1);
    allMeanSigW(ixN) = temp(2);
    allMeanSpW(ixN) = temp(3);
    
    temp2 = std(allParam(:,inx),0,2);
    allStdAveW(ixN) = temp2(1);
    allStdSigW(ixN) = temp2(2);
    allStdSpW(ixN) = temp2(3);
end
%}
dName = ['gmci_distri_sig_summData_',date,'.mat'];
save(fullfile(sFolder,dName),'allN','allM','allSp','allSig','allMeanFmin','allStdFmin',...
    'allMeanSpW','allStdSpW','allMeanAveW','allStdAveW','allMeanSigW','allStdSigW')


%% This part summarize data from Gcmi with correlated odor concentration, diff N
dataFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_corr';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
filesRaw = {allFile.name}';

% screen the files, only part are used. Since the folder may contain
% heterogeneous data set
str_marker = '2018-08-13';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

% define the basic parameters
allN = [20 30 40 50 100];
R = 10;      % number of receptors
num = 20;    % number of repeat
sp = 3;      % sparsity of input
sig = 2;  % different sigma

for i0 = 1:length(files)
    s1 = '(?<= *N)[\d.]+(?=_)';
%     s3 = '(?<= *_sig)[\d.]+(?=_)';
    N = str2num(char(regexp(files{i0},s1,'match')));
%     sig = str2num(char(regexp(files{i0},s3,'match')));
    inxN = find(allN == N);
%     inxSig = find(round(allSig*10) == round(sig*10));
    inxSig = 1;
    
    load(char(fullfile(dataFolder,filesep,files{i0})));
    if num ~= length(allfmin)
        disp('number of repeats in this simulation do not match!')
    end
    
    %hard threshold of sparisty
    thd = -4*sig;
    inx = allfmin < 2e2;  % only preserve the successed!
    allMeanFmin(inxSig,inxN) = mean(-allfmin(inx));
    allStdFmin(inxSig,inxN) = std(-allfmin(inx));
      
    meanW = [];
    stdW = [];
    for k0 = 1 :length(inx)
        temp = allMat(:,inx(k0));
        meanW = [meanW,mean(temp(temp > thd))];
        stdW = [stdW,std(temp(temp > thd))];     
    end
    allMeanAveW(inxSig,inxN) = mean(meanW);
    allMeanSigW(inxSig,inxN) = mean(stdW);
    allStdAveW(inxSig,inxN) = std(meanW);
    allStdSigW(inxSig,inxN) = std(stdW);
    
    temp2 = allMat(:,inx);
    allMeanSpW(inxSig,inxN) = mean(sum(temp2>thd,1)/size(allMat,1));
    allStdSpW(inxSig,inxN) = std(sum(temp2>thd,1)/size(allMat,1));

end
% save the data
saveName = ['gcmiCorr_diffN_M',num2str(R),'sp',num2str(sp),'sig',num2str(sig),'_',date,'.mat'];
dataName = fullfile(saveFolder,saveName);
save(dataName,'R','allN','sp','sig','allMeanFmin','allStdFmin','allMeanAveW',...
    'allStdAveW','allStdSigW','allStdSigW','allMeanSpW','allStdSpW')

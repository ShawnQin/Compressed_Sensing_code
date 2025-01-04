% analyze and plot the optimal interaction matrix
% this program load the folder that contain the simulation resulsts of
% opitmal W matrix, it plot the change of  target fucntion, the sparsity (fitted by a Gaussian)
% this is the second version of

close all
clear
clc

% define a globle variable, used in cmaes search

%% set the folder of data files, or use a gui to mannually specify the data location of data

dataFolder = '/Users/shan/Dropbox/olfactionProject/data/twoByMInt_diffSig_h2_05092018';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
filesRaw = {allFile.name}';

% screen the files, only part are used. Since the folder may contain
% heterogeneous data set
str_marker = '2018-05-09';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

% define the basic parameters
% R = [60,100,150,200];
num = 30;  % number of repeat
N = 2;     %number of odorants

% data structure store all the data
dataSumm = struct('allFmean',[],'allFstd',[],'R',[],'allSig',[],'spMean',[],'spStd',[]);

% string flag to extract information
Rstr = '(?<= *N2_R)[\d.]+(?=_)';
sigStr = '(?<= *_sig)[\d.]+(?=_)';


% cycle through all the data file, get R and allSig information
allR = [];
allSig = [];
for i0 = 1:length(files)
    allR = [allR,str2num(char(regexp(files{i0},Rstr,'match')))];
    allSig = [allSig, str2num(char(regexp(files{i0},sigStr,'match')))];
end
[val1,~,Rinx] = unique(allR);
[val2,~,Sinx] = unique(allSig);
dataSumm.R = sort(val1);
dataSumm.allSig = sort(val2);

% add a matrix of cell array to save all the matrix
dataSumm.allW = cell(length(dataSumm.R),length(dataSumm.allSig));

dataSumm.allFmean = zeros(length(dataSumm.R),length(dataSumm.allSig));
dataSumm.allFstd = zeros(length(dataSumm.R),length(dataSumm.allSig));
for i0 = 1:length(files)
    temp = load(char(fullfile(dataFolder,filesep,files{i0})));
    if num ~= length(temp.allfmin)
        error('number of repeats in this simulation do not match!')
    else
        sig = dataSumm.allSig(Sinx(i0));
        R = dataSumm.R(Rinx(i0));
        dataSumm.allFmean(Rinx(i0),Sinx(i0)) = mean(-temp.allfmin) + N/2*log(2*pi*exp(1)*sig^2);
        dataSumm.allFstd(Rinx(i0),Sinx(i0)) = std(temp.allfmin);
        
        %hard threshold of sparisty
        thd = -4*dataSumm.allSig(Sinx(i0));
%         tm = temp.allMat(:);
        sp = sum(temp.allMat > thd,1)/R/N;
%         dataSumm.sparsity(Rinx(i0),Sinx(i0)) = sum(tm > thd)/length(tm);
        dataSumm.spMean(Rinx(i0),Sinx(i0)) = mean(sp);
        dataSumm.spStd(Rinx(i0),Sinx(i0)) = std(sp);
        
        % save the matrix columwize
        m1 = temp.allMat(1:R,:);
        m2 = temp.allMat((R+1):end,:);
        dataSumm.allW{Rinx(i0),Sinx(i0)} = [m1(:),m2(:)];  %normalize

    end
 
end
% save the data
dataName = fullfile(saveFolder,'N2diffSigSummary_h2_05092018.mat');
save(dataName,'dataSumm')


%% plot the figures

%%%%%% IMPORTANT %%%%%%%%%
% set the Hill coefficient in response function 
h  = 2; 
%%%%%%%%%%%%%%%%%%%%%%%%


% set up the colors set that will be used
rdColor = brewermap(11,'RdBu');
buColor = brewermap(11,'Blues');

%%%%%% plot fmin as function of \sigma_c %%%%%%%%%%
figure
% hold on
% ah = cell(1,length(dataSumm.R));
errorbar(kron(dataSumm.allSig',ones(1,length(dataSumm.R))),dataSumm.allFmean',dataSumm.allFstd',...
    'o-','MarkerSize',12,'LineWidth',2)
% legend(cellstr(num2str(dataSumm.R)),'Location','northwest')
% legend('M = 60','M = 100','M = 150','M = 200','Location','northwest')
legend('M = 60','M = 100','M = 200','Location','northwest')

legend boxoff
% for i0 = 1:length(dataSumm.R)
%     ah = errorbar(dataSumm.allSig,dataSumm.allFmean(i0,:),dataSumm.allFstd(i0,:),...
%         'o-','MarkerSize',12,'LineWidth',2,'Color',buColor(i0*2 + 1,:));
% end
% hold off
% legend(ah,cellstr(num2str(dataSumm.R)),'Location','northwest')
% legend([fh1a fh1b],{'\sigma=2','\sigma=4'})

xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('diffeential entropy','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)

figPref = ['2xMDiffEntr_sig_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


%%%%%%%% plot the sparsity %%%%%%%%%%%%%%%
figure
colors = buColor(5:2:11,:);
set(groot,'defaultAxesColorOrder',colors)

errorbar(kron(dataSumm.allSig',ones(1,length(dataSumm.R))),dataSumm.spMean',dataSumm.spStd',...
    'o-','MarkerSize',12,'LineWidth',2)
% legend('M = 60','M = 100','M = 150','M = 200','Location','northwest')
legend('M = 60','M = 100','M = 200','Location','northwest')

legend boxoff

xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('sparsty of W','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)

figPref = ['2xMSparsity_sig_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


%% histogram and correlation of elements in W

% histogram of active elements as function of sigma
set(groot,'defaultAxesColorOrder','default')
figure

sInx = [1,4,7,11];  %select part to plot
for i0 = 1:length(dataSumm.R)
    subplot(2,2,i0)
%     for j0 = 1:10
    for j0 = 1:length(sInx)

        sig = dataSumm.allSig(sInx(j0));
        thd = -4*sig;
%         sig = dataSumm.allSig(j0);
        
%         X = -4*sig:0.02:4*sig;
        X = -4:0.02:4;        
        temp = dataSumm.allW{i0,sInx(j0)};
%         temp = dataSumm.allW{i0,j0};
   
        temp2 = temp(temp > thd)/sig;
%         histogram(temp2,'Normalization','probability')
        hold on
        [f,xi] = ksdensity(temp2,'Bandwidth',0.5); 
        plot(xi,f,'LineWidth',2);
        title(['M=',num2str(dataSumm.R(i0))])
%         pd = fitdist(temp2,'Kernel');
%         Y = pdf(pd,X);
%         plot(X,Y,'LineWidth',2)
        set(gca,'FontSize',18,'LineWidth',1,'XLim',[-4,4])
        xlabel('$\log(W/\sigma_c)$','Interpreter','latex','FontSize',16)
        ylabel('probability','FontSize',20)
%         hold off
    end
%     legend('\sigma = 1.5','\sigma = 1.9','\sigma = 2.5','\sigma = 3','Location','northwest')
    legend(['\sigma = ',num2str(allSig(sInx(1)))],['\sigma = ',num2str(allSig(sInx(2)))],...
        ['\sigma = ',num2str(allSig(sInx(3)))],['\sigma = ',num2str(allSig(sInx(4)))])
    legend boxoff
end

figPref = ['2xMHistoActi_sig_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

%%%%%%% plot the distribution of active element %%%%%%%%%%%%%
figure
for i0 = 1:length(dataSumm.R)
    subplot(2,2,i0)
    for j0 = 1:length(sInx)
        sig = dataSumm.allSig(sInx(j0));        
        temp = dataSumm.allW{i0,sInx(j0)};
        thd = -4*sig;
        inxA = temp>thd;
        inxB = temp>thd;
        
        temp2 = max(temp(xor(inxA(:,1),inxA(:,2)),:),[],2)/sig;
        hold on
        [f,xi] = ksdensity(temp2,'Bandwidth',0.3); 
        plot(xi,f,'LineWidth',2);
        
        set(gca,'FontSize',18,'LineWidth',1,'XLim',[-4,4])
        xlabel('$\log(W/\sigma_c)$','Interpreter','latex','FontSize',16)
        ylabel('probability','FontSize',20)
    end
    title(['M=',num2str(dataSumm.R(i0))])
    legend(['\sigma = ',num2str(allSig(sInx(1)))],['\sigma = ',num2str(allSig(sInx(2)))],...
        ['\sigma = ',num2str(allSig(sInx(3)))],['\sigma = ',num2str(allSig(sInx(4)))])
%     legend('\sigma = 1.5','\sigma = 1.9','\sigma = 2.5','\sigma = 3','Location','northwest')
    legend boxoff
end

figPref = ['2xMKernelActiSingle_sig_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

%%%%%%%% scatter plot of active elements %%%%%%%%
% sInx = [3,5,8,10];  %select part to plot
sInx = [1,4,7,11];

MAX_FIG = 100;
for i0 = 1:length(dataSumm.R)
    figure(MAX_FIG + i0)
    hold on
    for j0 = 1:length(sInx)
        sig = dataSumm.allSig(sInx(j0));
        temp = dataSumm.allW{i0,sInx(j0)}/sig;
        scatter(temp(:,1),temp(:,2),'LineWidth',1)  
        xlim([-4,-0.5])
        ylim([-4,-0.5])
    end
    
    set(gca,'FontSize',24,'LineWidth',1.5,'XLim',[-4,0.5],'YLim',[-4,0.5])
    title(['M=',num2str(dataSumm.R(i0))])
    xlabel('$w_1/\sigma_c$','Interpreter','latex','FontSize',28)
    ylabel('$w_2/\sigma_c$','Interpreter','latex','FontSize',28)
%     legend('\sigma = 1.5','\sigma = 1.9','\sigma = 2.5','\sigma = 3','Location','northeast')
    legend(['\sigma = ',num2str(allSig(sInx(1)))],['\sigma = ',num2str(allSig(sInx(2)))],...
        ['\sigma = ',num2str(allSig(sInx(3)))],['\sigma = ',num2str(allSig(sInx(4)))])
    legend boxoff
    hold off
    figPref = ['2xMScatterActi_sig_R',num2str(dataSumm.R(i0)),'_',date];
    saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
    print('-depsc',[saveFolder,filesep,figPref,'.eps'])
end


%%%%%%  analyze the ratio of both active part elements %%%%%%%%%
overLapMean = zeros(length(dataSumm.R),length(dataSumm.allSig));
overLapStd = zeros(length(dataSumm.R),length(dataSumm.allSig));
for i0 = 1:length(dataSumm.R)
    R = dataSumm.R(i0);
    for j0 = 1:length(dataSumm.allSig)
        temp = dataSumm.allW{i0,j0};
        overlp = zeros(num,1);
        thd = -4*dataSumm.allSig(j0);
        for k0 = 1:round(size(temp,1)/R)
            temp2 = temp((k0-1)*R+1:k0*R,:);
            inx = temp2(:,1) > thd & temp2(:,2) > thd;
            overlp(k0) = sum(inx)/R;
        end
        
        overLapMean(i0,j0) = mean(overlp);
        overLapStd(i0,j0) = std(overlp);
    end
    
end

% plot the overlap ratio
figure
colors = buColor(5:2:11,:);
set(groot,'defaultAxesColorOrder',colors)

errorbar(kron(dataSumm.allSig',ones(1,length(dataSumm.R))),overLapMean',overLapStd',...
    'o-','MarkerSize',12,'LineWidth',2)
hold on
% plot([1.1,3],[0.5,0.5],'--','LineWidth',2)
legend('M = 60','M = 100','M = 200','Location','southeast')
% legend('M = 60','M = 100','M = 150','M = 200','Location','southeast')

legend boxoff

xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('overlap ratio','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)

figPref = ['2xMActiOverlap_sig_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


%% scatter plot for larger sigma
set(groot,'defaultAxesColorOrder','default')
figure
hold on
for j0 = 1:3
    sig = dataSumm.allSig(j0);
    temp = dataSumm.allW{j0}/sig;
    scatter(temp(:,1),temp(:,2),'LineWidth',1)        
end
hold off
set(gca,'FontSize',24,'LineWidth',1.5,'XLim',[-4,0.5],'YLim',[-4,0.5])
title(['R=',num2str(dataSumm.R)])
xlabel('$w_1/\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('$w_2/\sigma_c$','Interpreter','latex','FontSize',28)
legend('\sigma = 3','\sigma = 4','\sigma = 5','Location','northeast')
legend boxoff

for i0 = 1:length(dataSumm.R)
    figure(MAX_FIG + i0)
    hold on
    for j0 = 1:length(sInx)
        sig = dataSumm.allSig(sInx(j0));
        temp = dataSumm.allW{i0,sInx(j0)}/sig;
        scatter(temp(:,1),temp(:,2),'LineWidth',1)        
    end
    set(gca,'FontSize',24,'LineWidth',1.5,'XLim',[-4,0.5],'YLim',[-4,0.5])
    title(['M=',num2str(dataSumm.R(i0))])
    xlabel('$w_1/\sigma_c$','Interpreter','latex','FontSize',28)
    ylabel('$w_2/\sigma_c$','Interpreter','latex','FontSize',28)
    legend(['\sigma = ',num2str(allSig(sInx(1)))],['\sigma = ',num2str(allSig(sInx(2)))],...
        ['\sigma = ',num2str(allSig(sInx(3)))],['\sigma = ',num2str(allSig(sInx(4)))])
%     legend('\sigma = 1.5','\sigma = 1.9','\sigma = 2.5','\sigma = 3','Location','northeast')
    legend boxoff
    hold off
    figPref = ['2xMScatterActi_sig_R',num2str(dataSumm.R(i0)),'_','h',num2str(h),'_',date];
    saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
    print('-depsc',[saveFolder,filesep,figPref,'.eps'])
end


%% fit a deterministic function between the w1 and w2
% first guess, it should have a form (|x|^p + |y|^p)^{1/p} = c0
% we can first guess p = 1/2
[ROW,COL] = size(dataSumm.allW);
dataSumm.p = cell(ROW,COL);
dataSumm.c0 = cell(ROW,COL);

for i0 = 3:ROW
    for j0 = 4:COL
        sig = dataSumm.allSig(j0);
        rawM = dataSumm.allW{i0,j0}/sig;
        [p,c0] = fitpNorm(rawM);
        dataSumm.p{i0,j0} = p;
        dataSumm.c0{i0,j0} = c0;
    end
end

% save the daa
dFile = fullfile(saveFolder,'summData_0509_Regul.mat');
% dFile = fullfile(saveFolder,'summDataSupp_0418_2.mat');
save(dFile,'-struct','dataSumm')


%% plot how the fitted parameter change
dFile = fullfile(saveFolder,'summData_0418_Regul.mat');
dataSumm = load(dFile);

% suplementary data
dataSupp = load(fullfile(saveFolder,'summDataSupp_0418_2.mat'));

% plot how p and c0 changes with sigma
ROW = 4;
COL = 10;
firstRingP = zeros(ROW,COL-1);  %the first column is negnegted
for i0 = 1:ROW
    for j0 = 2:COL
        firstRingP(i0,j0-1) = dataSumm.p{i0,j0}(1);
    end
end
figure
% meanP = cellfun(@mean,dataSumm.p(1));
plot(dataSumm.allSig(2:end),firstRingP,'o-','LineWidth',2,'MarkerSize',10)
ylim([0,1])
lg = legend('R=60','R=100','R=150','R=200','Location','northeast');
set(lg,'FontSize',20)
legend boxoff
xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('$p$','Interpreter','latex','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
figPref = ['2xMfitNormP_sig_Regul','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

%%%%%%% plot the c0 as a function of sigma_c %%%%%%
% the first ring
firstRing = zeros(ROW,COL-1);
for i0= 1:ROW
    for j0 = 2:COL
%         firstRing(i0,j0-1) = abs(dataSumm.c0{i0,j0}(1,1)).^(1/dataSumm.p{i0,j0}(1,1));
        firstRing(i0,j0-1) = dataSumm.c0{i0,j0}(1,1);
    end
end
plot(dataSumm.allSig(2:end),firstRing,'o-','LineWidth',2,'MarkerSize',10)
% ylim([0,1])
lg = legend('R=60','R=100','R=150','R=200','Location','northeast');
set(lg,'FontSize',20)
legend boxoff
title('First Ring')
xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('$c_0$','Interpreter','latex','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
figPref = ['2xMfitNormC0_sig_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


%%%% first ring, with data from supple, large sigma, only plot M = 100  %%
firstRingR100Sup = zeros(1,3);
for i0 = 1:3
    firstRingR100Sup(i0) = dataSupp.p{i0}(1);
end
firstRingR100 = [firstRingP(2,:),firstRingR100Sup];
sigR100 = [dataSumm.allSig(2:end),dataSupp.allSig];
figure
plot(sigR100,firstRingR100,'o-','LineWidth',3,'MarkerSize',12)
lg = legend('$R=100,|w_1|^p + |w_2|^p = c_0$','Location','northwest');
set(lg,'Interpreter','latex')
legend boxoff
xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('$c_0$','Interpreter','latex','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
figPref = ['2xMfitP_sig_R100_Regul','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

% the second ring, need data from the other data set
% first load the data from other data set
secondRing = zeros(1,5);
thirdRing = zeros(1,2);
selSig = [2.7,2.9,3,4,5];
for i0 = 1:2
    secondRing(i0) = dataSumm.c0{2,8+i0}(2);
end
for i0 = 1:3
    secondRing(2+i0) = dataSupp.c0{i0}(2);
end

for i0 = 1:2
    thirdRing(i0) = dataSupp.c0{1+i0}(3);
end
figure
hf1= plot(selSig,secondRing,'o-','LineWidth',2,'MarkerSize',12);
hold on
hf2 = plot(selSig(4:5),thirdRing,'^-','LineWidth',2,'MarkerSize',12);
hold off
lg = legend([hf1,hf2],'second ring','thrid ring','Location','southwest');
set(lg,'FontSize',20)
legend boxoff
xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('$c_0$','Interpreter','latex','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
figPref = ['2xMfitNormC0_scdThrdsig_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

% relative position of the three rings
figure
rings = [dataSupp.c0{2},dataSupp.c0{3}];
plot((1:3)',rings,'o-','LineWidth',2,'MarkerSize',12)
lg = legend('$\sigma_c = 4$','$\sigma_c = 5$','Location','northwest');
set(lg,'Interpreter','latex','FontSize',20)
legend boxoff
ylim([1.5,3])
xlabel('Ring index','FontSize',28)
ylabel('$c_0$','Interpreter','latex','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24,'XTick',1:1:3)
figPref = ['2xMfitNormC0_ringsPosi_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])



%%%%%%% plot the theoretical curves  %%%%%%%
allC0 = dataSupp.c0{2};
allP = dataSupp.p{2};
X = -3:0.02:0;
Y = zeros(length(allC0),length(X));
for i0 = 1:length(allC0)
    Y(i0,:) = -(allC0(i0) - abs(X).^allP(i0)).^(1/allP(i0));
end

figure
plot(X,Y,'LineWidth',3)
daspect([1 1 1])
ylim([-3,0])
xlabel('$w_1$','Interpreter','latex','FontSize',28)
ylabel('$w_2$','Interpreter','latex','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5)
figPref = ['2xMringsIllus_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


% plot one example of fit
p = dataSupp.p{1}(1);    % p
examplW = dataSupp.allW{1}/dataSupp.allSig(1);  % first ring
inx = abs(examplW(:,1)).^p + abs(examplW(:,2)).^p < 2.5;  % select data
plotData = examplW(inx,:);

figure
scatter(plotData(:,1),plotData(:,2),'k.')
hold on
c0 = dataSupp.c0{1}(1);
X = -3:0.01:0;
Y = -(c0 - abs(X).^p).^(1/p);
plot(X,Y,'r-','LineWidth',3)
xlim([-3,0])
ylim([-3,0])
xlabel('$w_1/\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('$w_2/\sigma_c$','Interpreter','latex','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
lgdStr = ['$|w_1|^{0.6}+ |w_1|^{0.6} =',num2str(c0),'$'];
lgd = legend('simulation',lgdStr,'Location','southwest');
set(lgd,'Interpreter','latex')
legend boxoff
figPref = ['2xMfitScatter_firstRing_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


%% density along the line
rawData = dataSupp.allW{2};
inx = all(rawData>-10,2);
scatter(rawData(inx,1),rawData(inx,2))

xlim([-10,0])
ylim([-10,0])

gridx1 = -10:.5:0;
gridx2 = -10:.5:0;
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
figure
ksdensity(rawData(inx,:),xi,'Bandwidth',0.4);
xlabel('$w_1$','Interpreter','latex','FontSize',28)
ylabel('$w_2$','Interpreter','latex','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5)

temp = rawData(inx,:);
pts = linspace(-10,0,25);
N = histcounts2(temp(:,1), temp(:,2), pts, pts);

% Create Gaussian filter matrix:
[xG, yG] = meshgrid(-12:1);
sigma = 1;
g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
g = g./sum(g(:));

figure
imagesc(pts, pts, conv2(N, g, 'same'));
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');



% project the data along the stretched curve
% select the data, only the first ring
p = dataSupp.p{1}(1);    % p
examplW = dataSupp.allW{1}/dataSupp.allSig(1);  % first ring
inx = abs(examplW(:,1)).^p + abs(examplW(:,2)).^p < 2.5;  % select data
plotData = examplW(inx,:);

% p = dataSupp.p{2}(1);
% inx = abs(plotData(:,1)).^p + abs(plotData(:,2)).^p < 2.5;  % select data
% X1 = rawData(inx,1);
% Y1 = rawData(inx,2);

X2 = plotData*[1;-1];
histogram(X2,'Normalization','probability')
xlabel('$|w|^p$','Interpreter','latex','FontSize',28)
ylabel('probability','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5)
figPref = ['2xM_wp_distAlongCurve','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

%% example
% Normally distributed sample points:
x = randn(1, 100);
y = randn(1, 100);

% Bin the data:
pts = linspace(-3, 3, 101);
N = histcounts2(y(:), x(:), pts, pts);

% Create Gaussian filter matrix:
[xG, yG] = meshgrid(-5:5);
sigma = 2.5;
g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
g = g./sum(g(:));

% Plot scattered data (for comparison):
subplot(1, 2, 1);
scatter(x, y, 'r.');
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));

% Plot heatmap:
subplot(1, 2, 2);
imagesc(pts, pts, conv2(N, g, 'same'));
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');



%%  ==================================================================== %%
%  This part analyze the optimal matrix simulated from gcmi method
%  For N = 100, M = 20, sig_c = 2 and SP = {1,2,3,4,5,8,10}
% ====================================================================

%% prepare and save the data, gcmi with different input sparsity (n)
dataFolder = '/Users/shan/Dropbox/olfactionProject/data/newGcmi20180831/GcmiSdp0828';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
filesRaw = {allFile.name}';

% screen the files, only part are used. Since the folder may contain
% heterogeneous data set
str_marker = '2018-08-29';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

% define the basic parameters
% R = [60,100,150,200];
R = 20;      % number of receptors
num = 40;    % number of repeat
N = 100;     % number of odorants
sig = 2;     % sigma_c
% h = 1;       % Hill coefficient of input function

% data structure store all the data
dataSumm = struct('allFmean',[],'allFstd',[],'sp',[],'allSig',[],'spMean',[],...
    'spStd',[],'meanAveW',[],'stdAveW',[],'meanStdW',[],'stdStdW',[]);
dataSumm.allSig = 2;  % sigma_c
% string flag to extract information
spStr = '(?<= *_S)[\d.]+(?=_)';

% cycle through all the data file, get R and allSig information
allSp = [];  %sparsity of input
for i0 = 1:length(files)
    allSp = [allSp, str2num(char(regexp(files{i0},spStr,'match')))];
end
% [val1,~,Rinx] = unique(allR);
[val2,~,Sinx] = unique(allSp);
dataSumm.sp = sort(val2);
% dataSumm.allSig = sort(val2);

% add a matrix of cell array to save all the matrix
dataSumm.allW = cell(length(dataSumm.sp),length(dataSumm.allSig));

dataSumm.allFmean = zeros(length(dataSumm.sp),length(dataSumm.allSig));
dataSumm.allFstd = zeros(length(dataSumm.sp),length(dataSumm.allSig));
for i0 = 1:length(files)
    temp = load(char(fullfile(dataFolder,filesep,files{i0})));
    if num ~= length(temp.allfmin)
        error('number of repeats in this simulation do not match!')
    else

        dataSumm.allFmean(Sinx(i0),1) = mean(-temp.allfmin) + N/2*log(2*pi*exp(1)*sig^2);
        dataSumm.allFstd(Sinx(i0),1) = std(temp.allfmin);
        
        %hard threshold of sparisty
        thd = -4*sig;
        
        sp = sum(temp.allMat > thd,1)/R/N;
        
        aveActW = zeros(size(temp.allMat,2),1);
        stdActW = zeros(size(temp.allMat,2),1);
        for k0= 1:size(temp.allMat,2)
%             aveActW(k0) = mean(temp.allMat(temp.allMat(:,k0) > thd,k0));
%             stdActW(k0) = std(temp.allMat(temp.allMat(:,k0) > thd,k0));
            t = temp.allMat(temp.allMat(:,k0) > thd,k0);
            actData = t(t<= -thd);
            pd = fitdist(actData,'normal');
            aveActW(k0) = pd.mean;
            stdActW(k0) = pd.sigma;
        end
        dataSumm.meanAveW(Sinx(i0),1) = mean(aveActW);
        dataSumm.stdAveW(Sinx(i0),1) = std(aveActW);
        dataSumm.meanStdW(Sinx(i0),1) = mean(stdActW);
        dataSumm.stdStdW(Sinx(i0),1) = std(stdActW);
%         dataSumm.sparsity(Rinx(i0),Sinx(i0)) = sum(tm > thd)/length(tm);
        dataSumm.spMean(Sinx(i0),1) = mean(sp);
        dataSumm.spStd(Sinx(i0),1) = std(sp);
        
        % save the matrix columwize
        dataSumm.allW{Sinx(i0),1} = temp.allMat;  %normalize

    end
 
end
% save the data
saveName = ['gcmi_Spdpd_N',num2str(N),'M',num2str(R),'Summary_',date,'.mat'];
dataName = fullfile(saveFolder,saveName);
save(dataName,'dataSumm')

%% plot the statistics of the matrix
% first, set the default graphics settings
defaultGraphicsSetttings

% the sparsity of W
figure
errorbar(dataSumm.sp',dataSumm.spMean,dataSumm.spStd,'o-','MarkerSize',12,'LineWidth',2)
legend('N=100,M=20,\sigma_c = 2, h = 1')
legend boxoff
xlabel('number of odors (sparsity)')
ylabel('sparsity of W')
figPref = ['N100M20_sparistyW_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


% mean and std of active elements
allMeanWmean = zeros(length(dataSumm.sp),1);
allMeanWstd = zeros(length(dataSumm.sp),1);
allStdWmean = zeros(length(dataSumm.sp),1);
allStdWstd = zeros(length(dataSumm.sp),1);

for i0 = 1:length(dataSumm.sp)
    thd = -5*sig;
    temp = dataSumm.allW{i0}(:);
    meanWeach = zeros(size(dataSumm.allW{i0},2),1);
    stdEach = zeros(size(dataSumm.allW{i0},2),1);
    for j0 = 1:size(dataSumm.allW{i0},2)
        temp2 = dataSumm.allW{i0}(:,j0);   
        meanWeach(j0) = mean(temp2(temp2 > thd));
        stdEach(j0) = std(temp2(temp2 > thd));
    end
    allMeanWmean(i0) = mean(temp(temp > thd));
    allMeanWstd(i0) = std(meanWeach);
    allStdWmean(i0) = mean(stdEach);
    allStdWstd(i0) = std(stdEach);
    
end

figure
errorbar(dataSumm.sp',allMeanWmean,allMeanWstd,'o-','MarkerSize',12,'LineWidth',2)
legend('N=100,M=20,\sigma_c = 2, h = 1')
legend boxoff
xlabel('number of odors (sparsity)')
ylabel('$\mu_w$','Interpreter','latex')
figPref = ['N100M20_meanActiveW_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

figure
errorbar(dataSumm.sp',allStdWmean,allStdWstd,'o-','MarkerSize',12,'LineWidth',2)
ylim([1.5,3])
legend('N=100,M=20,\sigma_c = 2, h = 1')
legend boxoff
xlabel('number of odors (sparsity)')
ylabel('$\sigma_w$','Interpreter','latex')
figPref = ['N100M20_sigW_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])



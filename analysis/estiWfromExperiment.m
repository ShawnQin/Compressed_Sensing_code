% this program estimate the sensitivity and plot the histogram of
% experimental W. The primiary data are from John Carlson and Guangwei's
% data
% The default Hill coefficient is 1, we can actually tune its value
% last revised in 07/16/2018
close all
clear

%% graphic settings, define color that might be used
defaultGraphicsSetttings
RdBu = brewermap(11,'RdBu');   % red and blue
Bu = brewermap(11,'Blues');    % blues
Gr = brewermap(11,'Greys');    % blues

lBu = [96,166,223]/255; %light blue
dpBu = [63,114,183]/255; % deep blue
dkBu = [50,78,147]/255;   %dark blue
Or = [220,150,71]/255;  % orange
brickRd = [201,69,89]/255;  %brick red
green = [107,169,81]/255;  %green
purple = [113,49,119]/255;  % purple

%% prepare the data and set up related functions
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
% Carlson's data
file1 = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/Carlson2006TableS1.xlsx';
[NUM1,~,~]=xlsread(file1);
file2 = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/CarslonORNdigit.xlsx';
[NUM2,~,~]=xlsread(file2);  %digital matrix defines strong excitation and inhibition


% for Carlson's data, we need to estimate the "true spontaneous activity"
% we define the "true" spontaneous activity as the maximum of two possible
% definition
r0 = max([NUM1(end,:);abs(min(NUM1,[],1))],[],1);  %adding 0.1 for stability
spMat = ones(110,1)*r0;
adjM = NUM1 + r0;  % adjust after considering the basal activity

% default hillCoef
h = 1;  % we can tune this to see how the reuslts change

% we assumed the same maxium spiking rate
Rmax = max(adjM(:)) + 1;  %to avoid numerical unstability, so we add 1

% alpha parameter
alp = Rmax./r0 - 1;  % each receptor has one alpha
alpMat = ones(110,1)*alp;  %this matrix is used for future use

%strong excitation based on the digital matrix
allM = adjM(1:end-1,:);
strongExi = allM(NUM2 > 0);
strongW = (alpMat(NUM2 > 0)./(Rmax./allM(NUM2 > 0) - 1) - 1).^(1/h);

% consider all excitation
allInx = allM > spMat;
allExci = allM(allInx);
exciW = (alpMat(allInx)./(Rmax./allM(allInx) - 1) - 1).^(1/h);

%strong inhibition
strongInhi = max(allM(NUM2 < 0),1);  % for stability
% strongInhi = allM(NUM2 < 0);  % for stability

strongInhiW = ((Rmax./strongInhi-1)./alpMat(NUM2 < 0) - 1).^(1/h);

% consider all excitation
allInx = allM < spMat;
allInhi = max(allM(allInx),1);   % for stability
% allInhi = allM(allInx);   % for stability

inhiW = ((Rmax./allInhi-1)./alpMat(allInx) - 1).^(1/h);

%% plot the figures
% ============================================================
% histogram of strong excitation and lognormal fit
% ============================================================
figure
hold on
set(gcf,'renderer','Painters')
histogram(log10(strongW),20,'Normalization','pdf')
xlim([-1,5])
% fit a lognormal distribution
pd = fitdist(log10(strongW),'normal');
X = -1:0.05:5;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y)
hold off
legend('excti w','Gassian fit')
legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
prefix = ['CarlsonStrongExcitW_fit_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% ============================================================
% histogram of all excitation and lognormal fit
% ============================================================
figure
hold on
set(gcf,'renderer','Painters')
histogram(log10(exciW),20,'Normalization','pdf')
xlim([-2,4])
% fit a lognormal distribution
pd = fitdist(log10(exciW),'normal');
X = -2:0.05:4;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y)
hold off
legend('excti w','Gassian fit')
legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
prefix = ['CarlsonAllExcitW_fit_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ============================================================
% histogram of strong inhibition and lognormal fit
% ============================================================
figure
hold on
set(gcf,'renderer','Painters')
histogram(log10(strongInhiW),20,'Normalization','pdf')
xlim([-1,2])
% fit a lognormal distribution
pd = fitdist(log10(strongInhiW),'normal');
X = -1:0.05:2;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y)
hold off
legend('inhi w','Gassian fit')
legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
prefix = ['CarlsonStrongInhiW_fit_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ============================================================
% histogram of all inhibition and lognormal fit
% ============================================================
figure
hold on
set(gcf,'renderer','Painters')
histogram(log10(inhiW),20,'Normalization','pdf')
xlim([-2,2])
% fit a lognormal distribution
pd = fitdist(log10(inhiW),'normal');
X = -2:0.05:2;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y)
hold off
legend('inhi w','Gassian fit')
legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
prefix = ['CarlsonAllInhiW_fit_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ============================================================
% merge histogram of all excitation and inhibition and lognormal fit
% ============================================================
figure
hold on
set(gcf,'renderer','Painters')
histogram(log10(exciW),20,'Normalization','count')
histogram(log10(inhiW),20,'Normalization','count')
hold off
legend('exci','inhi')
legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('frequency')
prefix = ['CarlsonMergeExciInhi_histo_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% ============================================================
% merged tight histogram of all excitation and inhibition and lognormal fit
% ============================================================
% histogram
figureSize = [0 0 4.4 4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
[ha, pos] = tight_subplot(2,1,[0 0],[.2 .01],[.2 .03]);

axes(ha(1))
h1 = histfit(log(exciW),20);
xlim([-5,8])
ylim([0,350])
set(gca,'XTick',[],'YTick',0:150:300,'FontSize',22);

axes(ha(2))
h2 =  histfit(log(inhiW),20);
xlim([-5,8])
ylim([0,350])
xlabel('$\ln(w)$','Interpreter','latex')
set(gca,'XTick',-5:5:5,'YTick',0:150:300,'FontSize',22);

% set the color
h1(1).FaceColor = Bu(5,:);h2(1).FaceColor = RdBu(4,:);
h1(2).Color = Bu(10,:);h2(2).Color =RdBu(2,:);
h1(2).LineWidth = 3;h2(2).LineWidth = 3;

figNamePref = ['Figure6_ExciInhi_expermient_histo_',date];
saveas(gcf,[saveFolder,filesep,figNamePref,'.fig'])
print('-depsc',[saveFolder,filesep,figNamePref,'.eps'])


%% Guangwei's old data
file3 = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/Guangwei_Log_10EC_50.xlsx';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

[NUM3,TEXT,~]=xlsread(file3);
allW = NUM3(abs(NUM3) > 0);

% ================================================
% plot the heatmap matrix
% ================================================
newM = -NUM3;
newM(:,[6,13,15]) = [];
odorLabel = cell(1,18);
receptorLabel = cell(1,18);
inx1 = [2:5,7:12,14,16:size(TEXT,2)];
for i0 = 1:length(inx1)
    receptorLabel{i0} = TEXT{1,i0};
end

for i0 = 1:18
    odorLabel{i0} = TEXT{i0+1,1};
end

figure;  pos = get(gcf, 'pos'); set(gcf, 'pos', [pos(1), pos(2), 610, 300]);
imagesc(newM); 
set(gca, 'CLim', [0 max(newM(:))]);
set(gca,'XTick',1:size(newM,2),'FontSize',16);
set(gca,'XTickLabel',receptorLabel);
set(gca,'xaxisLocation','top');
set(gca,'YTick',1:size(newM,1),'FontSize',16);
set(gca,'YTickLabel',odorLabel);
set(gca, 'XTickLabelRotation', 45);
colormap(jet);
c = colorbar; 
c.TickLabels{1} = 'NaN'; 
c.Label.String = '-log10(EC50)';
c.FontSize = 16;
prefix = ['Guangwei_EC50_heatMap_disorder',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ================================
% statistical test
% ================================
[H, pValue, SWstatistic] = swtest(newM(newM>0), 0.01);


% ================================================
% fit histogram and fit with Gaussian distribution
% ================================================
% 
figure
figureSize = [0 0 5 4.5];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
hold on
set(gcf,'renderer','Painters')
hh = histogram(allW,15,'Normalization','pdf',...
    'EdgeColor','none','FaceColor',lBu,'FaceAlpha',0.4);
stairs([hh.BinEdges(1),hh.BinEdges,hh.BinEdges(end)],...
    [0,hh.Values,hh.Values(end),0],'Color',dkBu,'LineWidth',2)

pd = fitdist(allW,'normal');
X = -8:0.05:-2;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y,'LineWidth',4,'Color',Or)
box on
hold off
% legend('w','Gassian fit')
% legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
set(gca,'Layer', 'top')
prefix = ['GuangweiStrongW_fit_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% check if it is random, this can be compared with our simulation


% ===============================================================
% Histogram comparison of the original and shuffled
% ===============================================================

w0 = NUM3;
w0(:,[6 13 15])= [];
corrType = 'Spearman';
L = size(w0,1)*size(w0,2);

rowCol = 'column';  % row or column

% original
if strcmp(rowCol,'row')
    w = w0';
elseif strcmp(rowCol,'column')
    w = w0;
end
    
t1 = corr(w,'type',corrType);
y1 = triu(t1,1);
[bc1, edg1] = histcounts(y1(y1~=0),15);

% shuffled 
newW = reshape(w(randperm(L)),size(w,1),size(w,2));
if strcmp(rowCol,'row')
    t2 = corr(newW','type',corrType);
elseif strcmp(rowCol,'column')
    t2 = corr(newW,'type',corrType);
end

y2 = triu(t2,1);
% [bc2, bp2] = hist(y2(y2~=0),15);
[bc2,edg2] = histcounts(y2(y2~=0),15);


figure
hold on
Y1 = [bc1;bc1]*2/size(w,1)/(size(w,1)-1); %normalized, probability
ar1 = area(sort(edg1([1:end-1 2:end])),Y1(1:end));
ar1.FaceAlpha = 0.3;
ar1.FaceColor = Or;
ar1.EdgeColor = Or;
ar1.LineWidth = 2;

Y2 = [bc2;bc2]*2/size(w,1)/(size(w,1)-1); %normalized, probability
ar2 = area(sort(edg2([1:end-1 2:end])),Y2(1:end));
ar2.FaceAlpha = 0.4;
ar2.FaceColor = Gr(4,:);
ar2.EdgeColor = Gr(8,:);
ar2.LineWidth = 2;
set(gca,'Layer','top')
hold off
xlim([-1,1])
lg = legend('original','shuffled');
set(lg,'FontSize',16)
legend boxoff
box on
xlabel('$\rho_s$','Interpreter','latex')
ylabel('probability')
xlim([-1,1])
set(gca,'XTick',-1:0.5:1)
prefix = ['Guangwei_',rowCol,'Shuff_',corrType,'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-painters','-dpdf',[saveFolder,filesep,prefix,'.pdf'])


% t1 = corr(w,'type',corrType);
% y1 = triu(t1,1);
rowMean = mean(y1(y1~=0));
rowStd = std(y1(y1~=0));

repeats = 100;

L = size(w,1)*size(w,2);
allMeanCorr = zeros(repeats,1);
allStdCorr = zeros(repeats,1);

for i0 = 1:repeats
    newW = reshape(w0(randperm(L)),size(w0,1),size(w0,2));
    if strcmp(rowCol,'row')
        temp = corr(newW','type',corrType);
    elseif strcmp(rowCol,'column')
        temp = corr(newW,'type',corrType);
    end
    y3 = triu(temp,1);
    allMeanCorr(i0) = mean(y3(y3~=0));
    allStdCorr(i0) = std(y3(y3~=0));
end

figure
hold on
box on
scatter(allMeanCorr,allStdCorr,20,Gr(8,:),'filled')
% xlim([-0.075,0.075])

scatter(rowMean,rowStd,80,Or,'filled')
hold off
lg = legend('shuffled','original');
xlabel('$\langle \rho_s \rangle$','Interpreter','latex')
ylabel('$\sigma_{\rho_s}$','Interpreter','latex')
prefix = ['Guangwei_',rowCol,'ShuffScatter_',corrType,'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% ========================================
% We need to show the skewness and Kurtosis of this data set
% ========================================
mu = mean(allW);
sd = std(allW);
skLarva = skewness(allW);
ktLarva = kurtosis(allW);

% compare with Gaussian distribution
NS = 200; %sample 200 gaussian distribtion with same mu and std
GaussSk = zeros(NS,1);
GaussKt = zeros(NS,1);
for i0 = 1:NS
    rd = randn(length(allW),1)*sd + mu;
    GaussSk(i0) = skewness(rd);
    GaussKt(i0) = kurtosis(rd);
end

%% Guangwei's new data set, paper acepted by Neuron

file3 = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/flyLarva_Guangwei_LogEC50_New.xlsx';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

[NUM3,TEXT,~]=xlsread(file3);
allW = NUM3(abs(NUM3) > 0);

% ================================================
% plot the heatmap matrix
% ================================================
newM = -NUM3;
% newM(:,[6,13,15]) = [];
% odorLabel = cell(1,18);
% receptorLable = cell(1,18);
% inx1 = [2:5,7:12,14,16:size(TEXT,2)];
% for i0 = 1:length(inx1)
%     receptorLable{i0} = TEXT{1,i0};
% end
receptorLabel = TEXT(1,2:end);
odorLabel = TEXT(2:end,1);

sr = '[\w-\s,]{3,30}';
for i0 = 1:length(receptorLabel)
    receptorLabel{i0} = char(regexp(receptorLabel{i0},sr,'match'));
end

for i0 = 1:length(odorLabel)
    odorLabel{i0} = char(regexp(odorLabel{i0},sr,'match'));
end
% for i0 = 1:18
%     odorLabel{i0} = TEXT{i0+1,1};
% end

figure;  pos = get(gcf, 'pos'); set(gcf, 'pos', [pos(1), pos(2), 8, 6]);
imagesc(newM); 
set(gca, 'CLim', [0 max(newM(:))]);
set(gca,'XTick',1:size(newM,2),'FontSize',16);
set(gca,'XTickLabel',receptorLabel);
set(gca,'xaxisLocation','top');
set(gca,'YTick',1:size(newM,1),'FontSize',16);
set(gca,'YTickLabel',odorLabel);
set(gca, 'XTickLabelRotation', 45);
colormap(jet);
c = colorbar; 
c.TickLabels{1} = 'NaN'; 
c.Label.String = '-log10(EC50)';
c.FontSize = 16;
prefix = ['Guangwei_EC50_heatMap_disorder',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ================================
% statistical test
% ================================
[H, pValue, SWstatistic] = swtest(newM(newM>0), 0.01);


% ===============================================================
% Fit the data with a skewed gaussian
% ===============================================================
gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
skewStand = @(x,alpha) 2*gaussian(x).*normcdf(alpha*x);
skewPdf =  @(x,xi,omi,alp) 2/omi*gaussian((x-xi)/omi).*normcdf(alp*(x-xi)/omi);

dVec = newM(~isnan(newM));
% ecdf
[Y,X] = ecdf(dVec);

fun = @(x,xdata) normcdf((xdata-x(1))/x(2)) - 2*myOwenT((xdata-x(1))/x(2),x(3));
x0 = [mean(dVec),std(dVec),3];
lb = [0,1,0];
ub = [10,10,10];
% least square fit
optParam = lsqcurvefit(fun,x0,X(2:end),Y(2:end),lb,ub);

xi = optParam(1); omi = optParam(2); alp = optParam(3);

% ================================================
% fit histogram and fit with Gaussian distribution
% ================================================
figure
figureSize = [0 0 5 4.5];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
hold on
set(gcf,'renderer','Painters')
hh = histogram(abs(-allW),15,'Normalization','pdf',...
    'EdgeColor','none','FaceColor',lBu,'FaceAlpha',0.4);
% stairs([hh.BinEdges(1),hh.BinEdges,hh.BinEdges(end)],...
%     [0,hh.Values,hh.Values(end),0],'Color',dkBu,'LineWidth',2)

pd = fitdist(abs(allW),'normal');
X = 0:0.05:10;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y,'LineWidth',3,'Color',Or)
Y2 = skewPdf(X,xi,omi,alp);
plot(X,Y2,'LineWidth',3,'Color',RdBu(2,:))
box on
ylim([0,0.5])
hold off
lg = legend('experiment','Gaussian fit','Skewed Gaussian fit');
set(lg,'FontSize',16)
% legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
set(gca,'Layer', 'top')
prefix = ['GuangweiStrongW_fit_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])



% plot and compare
figure
hold on
plot(X,Y,'Color',Gr(6,:))
skewCDF = @(x) normcdf((x-xi)/omi) - 2*myOwenT((x-xi)/omi,alp);
plot(X(2:end),skewCDF(X(2:end)),'Color',Bu(9,:))
hold off
box on
lg = legend('experiment','Skewed Gaussian')
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('CDF')


% ===============================================================
% Histogram comparison of the original and shuffled
% ===============================================================
w0 = -NUM3;
% w0(:,[6 13 15])= [];
w0(isnan(w0)) = min(0,min(w0(:)));
corrType = 'Kendall';
L = size(w0,1)*size(w0,2);

rowCol = 'column';  % row or column

% original
if strcmp(rowCol,'row')
    w = w0';
elseif strcmp(rowCol,'column')
    w = w0;
end
    
t1 = corr(w,'type',corrType);
y1 = triu(t1,1);
[bc1, edg1] = histcounts(y1(y1~=0),15);

% shuffled 
newW = reshape(w(randperm(L)),size(w,1),size(w,2));
if strcmp(rowCol,'row')
    t2 = corr(newW','type',corrType);
elseif strcmp(rowCol,'column')
    t2 = corr(newW,'type',corrType);
end

y2 = triu(t2,1);
% [bc2, bp2] = hist(y2(y2~=0),15);
[bc2,edg2] = histcounts(y2(y2~=0),15);


figure
hold on
% Y1 = [bc1;bc1]*2/size(w,1)/(size(w,1)-1); %normalized, probability
Y1 = [bc1;bc1]/sum(bc1)/(edg1(2)- edg1(1)); %normalized, probability
ar1 = area(sort(edg1([1:end-1 2:end])),Y1(1:end));
ar1.FaceAlpha = 0.3;
ar1.FaceColor = Or;
ar1.EdgeColor = Or;
ar1.LineWidth = 2;

% Y2 = [bc2;bc2]*2/size(w,1)/(size(w,1)-1); %normalized, probability
Y2 =  [bc2;bc2]/sum(bc2)/(edg2(2)- edg2(1));
ar2 = area(sort(edg2([1:end-1 2:end])),Y2(1:end));
ar2.FaceAlpha = 0.4;
ar2.FaceColor = Gr(4,:);
ar2.EdgeColor = Gr(8,:);
ar2.LineWidth = 2;
set(gca,'Layer','top')
hold off
xlim([-1,1])
lg = legend('original','shuffled');
set(lg,'FontSize',16)
legend boxoff
box on
xlabel('$\rho_s$','Interpreter','latex')
ylabel('probability')
xlim([-1,1])
set(gca,'XTick',-1:0.5:1)
prefix = ['Guangwei_New',rowCol,'Shuff_',corrType,'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-painters','-dpdf',[saveFolder,filesep,prefix,'.pdf'])


% t1 = corr(w,'type',corrType);
% y1 = triu(t1,1);
rowMean = mean(y1(y1~=0));
rowStd = std(y1(y1~=0));

repeats = 100;

L = size(w,1)*size(w,2);
allMeanCorr = zeros(repeats,1);
allStdCorr = zeros(repeats,1);

for i0 = 1:repeats
    newW = reshape(w0(randperm(L)),size(w0,1),size(w0,2));
    if strcmp(rowCol,'row')
        temp = corr(newW','type',corrType);
    elseif strcmp(rowCol,'column')
        temp = corr(newW,'type',corrType);
    end
    y3 = triu(temp,1);
    allMeanCorr(i0) = mean(y3(y3~=0));
    allStdCorr(i0) = std(y3(y3~=0));
end

figure
hold on
box on
scatter(allMeanCorr,allStdCorr,20,Gr(8,:),'filled')
% xlim([-0.075,0.075])

scatter(rowMean,rowStd,80,Or,'filled')
hold off
lg = legend('shuffled','original');
xlabel('$\langle \rho_s \rangle$','Interpreter','latex')
ylabel('$\sigma_{\rho_s}$','Interpreter','latex')
prefix = ['Guangwei_New_',rowCol,'ShuffScatter_',corrType,'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% ========================================
% We need to show the skewness and Kurtosis of this data set
% ========================================
mu = mean(-allW);
sd = std(allW);
skLarva = skewness(-allW);
ktLarva = kurtosis(-allW);

% compare with Gaussian distribution
NS = 1000; %sample 200 gaussian distribtion with same mu and std
GaussSk = zeros(NS,1);
GaussKt = zeros(NS,1);
for i0 = 1:NS
    rd = randn(length(allW),1)*sd + mu;
    GaussSk(i0) = skewness(rd);
    GaussKt(i0) = kurtosis(rd);
end

% make a histogram to compare the Gaussian distribution and the experiment
% first, compare the sknewness
figure
% set(gcf,'Renderer','painters')
hold on
histogram(GaussSk,'FaceColor',Bu(9,:),'Normalization','probability');
ah = gca;
ylim = ah.YLim;
plot([skLarva; skLarva],ylim,'r--','LineWidth',2)
box on
set(gca,'XLim',[-0.5, 1])
legend('Gaussian','Experiment')
xlabel('skewness')
ylabel('probability')
prefix = ['Guangwei_Skewness_gauss_',corrType,'_',trc,'_N',num2str(N),'M',...
    num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


%second, compare the excess kurtosis, defined as original minus 3
figure
% set(gcf,'Renderer','painters')
hold on
histogram(GaussKt-3,'FaceColor',Bu(9,:),'Normalization','probability');
ah = gca;
ylim = ah.YLim;
plot([ktLarva-3; ktLarva-3],ylim,'r--','LineWidth',2)
box on
set(gca,'XLim',[-1, 1.6])
legend('Gaussian','Experiment')
xlabel('excess kurtosis')
ylabel('probability')
prefix = ['Guangwei_Kutosis_gauss_',corrType,'_',trc,'_N',num2str(N),'M',...
    num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


%% Saito 2009 data
% in this data set contain a panel of 63 odorants, 53 mouse Or and 10 Human
% Ors
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
dataFile = "/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/Saito2009s3.xls";
[NUM4,TEXT,~]=xlsread(dataFile);

% humman receptor index
humanInx = [1,5,14,17,18,20,23,34,35,40];
allInx = 1:1:62;
selInx = setdiff(allInx,humanInx);
mouseMat = NUM4(selInx,:);  % mouse-only data
mouseVec = mouseMat(mouseMat ~= 0);
allVec = NUM4(NUM4 ~=0);

% name of mouse receptor
receptorLabel = cell(length(selInx),1);
odorLabel = cell(size(mouseMat,2),1);
for i0 = 1:length(selInx)
    str = TEXT{1+selInx(i0),1};
    out = regexp(str, '''', 'split');
    receptorLabel{i0} = char(out(2:2:end));
end

% name of odorants
for i0 = 1:size(mouseMat,2)
    str = TEXT{1,1+i0};
    out = regexp(str, '''', 'split');
    odorLabel{i0} = char(out(2:2:end));
end
% ====================================
% heat map of EC50
% ====================================
figure;  pos = get(gcf, 'pos'); set(gcf, 'pos', [pos(1), pos(2), 8, 12]);
imagesc(abs(mouseMat')); 
set(gca, 'CLim', [0 max(abs(mouseMat(:)))]);
set(gca,'XTick',1:size(mouseMat,1),'FontSize',12);
set(gca,'XTickLabel',receptorLabel);
set(gca,'xaxisLocation','top');
set(gca,'YTick',1:size(mouseMat,2),'FontSize',12);
set(gca,'YTickLabel',odorLabel);
set(gca, 'XTickLabelRotation', 45);
colormap(jet);
c = colorbar; 
c.TickLabels{1} = 'NaN'; 
c.Label.String = '-log10(EC50)';
c.FontSize = 12;

prefix = ['Saito_Mouse_EC50_heatMap_disorder',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])
% ================================
% statistical test, test if data is Gaussian
% ================================
[H, pValue, SWstatistic] = swtest(allVec, 0.01);


% plot the histogram of mouse only data and add a Gaussian fit
figure
figureSize = [0 0 5 4.5];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
hold on
set(gcf,'renderer','Painters')
h1 = histogram(-mouseVec,25,'Normalization','pdf','FaceColor',lBu,'FaceAlpha',...
    0.4,'EdgeColor','none');
stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
    [0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',2)
% xlim([-2,2])
% fit a lognormal distribution
pd = fitdist(-mouseVec,'normal');
X = 2:0.05:8;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y,'LineWidth',4,'Color',Or)
hold off
lg = legend('Experiment','Normal distribution','Location','northwest');
set(lg,'FontSize',16)
legend boxoff
box on

% legend('w','Gassian fit')
% legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
prefix = ['Saito_mouse_fit_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ================================================
% plot all the EC50, both mouse and human data are used
% ================================================
figure
hold on
set(gcf,'renderer','Painters')
histogram(allVec,25,'Normalization','pdf')
% xlim([-2,2])
% fit a lognormal distribution
pd = fitdist(allVec,'normal');
X = -7:0.05:-2;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y,'LineWidth',4,'Color',RdBu(3,:))
hold off
lg = legend('Experiment','Normal distribution','Location','northwest');
set(lg,'FontSize',16)
legend boxoff
box on

% legend('w','Gassian fit')
% legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
prefix = ['Saito_allData_fit_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


%% Mainland 2015 data
% these two data set is extracted from R program
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
dataFile1 = "/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/Mainland2015/	.xls";
dataFile2 = "/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/Mainland2015/EC50StrigentCon_Mainland2015.xls";

[NUM5,~,~]=xlsread(dataFile1);
[NUM6,~,~]=xlsread(dataFile2);


% sensitivity in a loose criterior,
Y1 = 1./NUM5(~isnan(NUM5));

% sensitivity in a strigent criterior
Y2 = 1./NUM6;

% ===============================================
% Statistical test,Jarque-Bera test
% ==========================================
% [flg1,p1,jbstat1,critval1] = jbtest(log10(Y1),[],0.0001);
[H1, pValue1, SWstatistic1] = swtest(log10(Y1), 0.01);
% [flg2,p2,jbstat2,critval2] = jbtest(log10(Y2),[],0.0001);
[H2, pValue2, SWstatistic2] = swtest(log10(Y2), 0.01);
% ==========================================
% histogram of sensitivity matrix, loose criterion
% ==========================================
figure
hold on
set(gcf,'renderer','Painters')
% histogram(log10(Y1),20,'Normalization','pdf')
h1 = histogram(log10(Y1),20,'Normalization','pdf','FaceColor',lBu,'FaceAlpha',...
    0.4,'EdgeColor','none');
stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
    [0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',2)
% xlim([-2,2])
% fit a lognormal distribution
pd = fitdist(log10(Y1),'normal');
X = 0:0.1:10;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y,'LineWidth',4,'Color',Or)
hold off
lg = legend('Experiment','Normal distribution','Location','northeast');
set(lg,'FontSize',16)
legend boxoff
box on

% legend('w','Gassian fit')
% legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
prefix = ['Mainland_human_loose_fit_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ==========================================
% histogram of sensitivity matrix, stringent criterion
% ==========================================
figure
hold on
set(gcf,'renderer','Painters')
histogram(log10(Y2),12,'Normalization','pdf')
% xlim([-2,2])
% fit a lognormal distribution
pd = fitdist(log10(Y2),'normal');
X = 0:0.1:10;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y,'LineWidth',4,'Color',RdBu(3,:))
hold off
lg = legend('Experiment','Normal distribution','Location','northeast');
set(lg,'FontSize',16)
legend boxoff
box on

% legend('w','Gassian fit')
% legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
prefix = ['Mainland_human_stringent_fit_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


%% Carlson's 2008 Neuron paper, fly larvae
% 21 Ors and 27 odorants, with basal activity
% concentration: 10^-2 dilution

saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
dataFile = "/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/Carlson2008.xlsx";
[NUM5,TEXT,~] = xlsread(dataFile);
nOR = 21;
nOdor = 27;

% change of firing rate matrix
MatFRchg = NUM5(1:nOdor,:);

% spontaneus activity
spon = NUM5(nOdor+1,:);

% for Carlson's data, we need to estimate the "true spontaneous activity"
% we define the "true" spontaneous activity as the maximum of two possible
% definition
r0 = max([spon;abs(min(MatFRchg,[],1))],[],1);  %adding 0.1 for stability
spMat = ones(nOdor,1)*r0;
allM = MatFRchg + r0;  % adjust after considering the basal activity

% default hillCoef
h = 1;  % we can tune this to see how the reuslts change

% we assumed the same maxium spiking rate
Rmax = max(allM(:)) + 1;  %to avoid numerical unstability, so we add 1

% alpha parameter
alp = Rmax./r0 - 1;  % each receptor has one alpha
alpMat = ones(nOdor,1)*alp;  %this matrix is used for future use

%strong excitation based on the digital matrix
% allM = adjM(1:end-1,:);
% strongExi = allM(NUM2 > 0);
% strongW = (alpMat(NUM2 > 0)./(Rmax./allM(NUM2 > 0) - 1) - 1).^(1/h);

% consider all excitation
allInx = allM > spMat;
allExci = allM(allInx);
exciW = (alpMat(allInx)./(Rmax./allM(allInx) - 1) - 1).^(1/h);

%strong inhibition
% strongInhi = max(allM(NUM2 < 0),1);  % for stability
% % strongInhi = allM(NUM2 < 0);  % for stability
% strongInhiW = ((Rmax./strongInhi-1)./alpMat(NUM2 < 0) - 1).^(1/h);

% consider all inhibition
allInx = allM < spMat;
allInhi = max(allM(allInx),1);   % for stability
% allInhi = allM(allInx);   % for stability

inhiW = ((Rmax./allInhi-1)./alpMat(allInx) - 1).^(1/h);

% plot the histogram
% histogram
figureSize = [0 0 4.5 4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
[ha, pos] = tight_subplot(2,1,[0 0],[.2 .03],[.2 .03]);

axes(ha(1))
h1 = histfit(log(exciW),20);
xlim([-5,8])
ylim([0,60])
set(gca,'XTick',[],'YTick',0:20:60,'FontSize',22);

axes(ha(2))
h2 =  histfit(log(inhiW),20);
xlim([-5,8])
ylim([0,60])
xlabel('$\ln(w)$','Interpreter','latex')
set(gca,'XTick',-5:5:5,'YTick',0:20:60,'FontSize',22);

% set the color
h1(1).FaceColor = Bu(5,:);h2(1).FaceColor = RdBu(4,:);
h1(2).Color = Bu(10,:);h2(2).Color =RdBu(2,:);
h1(2).LineWidth = 3;h2(2).LineWidth = 3;

figNamePref = ['Carlson2008Larvae_ExciInhi_histo_',date];
saveas(gcf,[saveFolder,filesep,figNamePref,'.fig'])
print('-depsc',[saveFolder,filesep,figNamePref,'.eps'])


% ======================================================================
% calculate the fraction of strong excitation and inhibtion with new
% criterior, excitation: > 2 fold spontaneous rate; inihibition: < 50% of
% basal activity
% ======================================================================
spMat = ones(nOdor,1)*spon;
allM = MatFRchg + spon;  % adjust after considering the basal activity
allM(allM<0) = 0;

Rmax = max(allM(:)) + 1;  %to avoid numerical unstability, so we add 1
alp = Rmax./spon - 1;  % each receptor has one alpha
alpMat = ones(nOdor,1)*alp;  %this matrix is used for future use

allInx2 = allM > 2*spMat;
allExci2 = allM(allInx2);
exciW2 = (alpMat(allInx2)./(Rmax./allM(allInx2) - 1) - 1).^(1/h);

allInx2 = allM < spMat/2;
allInhi2 = max(allM(allInx2),0.1);   % for stability
% allInhi = allM(allInx);   % for stability


inhiW2 = ((Rmax./allInhi2-1)./alpMat(allInx2) - 1).^(1/h);


% plot a merged histogram
figureSize = [0 0 4.5 4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
[ha, pos] = tight_subplot(2,1,[0 0],[.2 .03],[.2 .03]);

axes(ha(1))
h1 = histfit(log(exciW2),20);
xlim([-5,8])
ylim([0,60])
set(gca,'XTick',[],'YTick',0:20:60,'FontSize',22);

axes(ha(2))
h2 =  histfit(log(inhiW2),20);
xlim([-5,8])
ylim([0,60])
xlabel('$\ln(w)$','Interpreter','latex')
set(gca,'XTick',-5:5:5,'YTick',0:20:60,'FontSize',22);

% set the color
h1(1).FaceColor = Bu(5,:);h2(1).FaceColor = RdBu(4,:);
h1(2).Color = Bu(10,:);h2(2).Color =RdBu(2,:);
h1(2).LineWidth = 3;h2(2).LineWidth = 3;

figNamePref = ['Carlson2008Larvae_ExciInhi_histo_',date];
saveas(gcf,[saveFolder,filesep,figNamePref,'.fig'])
print('-depsc',[saveFolder,filesep,figNamePref,'.eps'])


%% Compare Calcium imaging and electrophysiology, Guangwei's and Kreher's data

%load the data set
file1 = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/flyLarva_Guangwei_LogEC50_New.xlsx';
file2 = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/Carlson2008.xlsx';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% load Guangwei's data
[NUM1,TEXT1,~]=xlsread(file1);
allW1 = NUM1(abs(NUM1) > 0);

% load Carlson's data
[NUM2,TEXT2,~]=xlsread(file2,1);   %sheet one is the data with highest concentration
[NUM4,TEXT4,~]=xlsread(file2,2);   %sheet one is the data with highest concentration
[sem1,~,~]=xlsread(file2,2);   %sheet one is the data with highest concentration
[sem2,~,~]=xlsread(file2,2);   %sheet one is the data with highest concentration

NUM2(27,:) = [];  % delete the CO2 data, since it only apper in this matrix
allW2 = NUM2(abs(NUM2) > 0);

receptorLabel = TEXT1(1,2:end);
odorLabel = TEXT1(2:end,1);

sr = '[\w-\s,]{3,30}';
for i0 = 1:length(receptorLabel)
    receptorLabel{i0} = char(regexp(receptorLabel{i0},sr,'match'));
end

for i0 = 1:length(odorLabel)
    odorLabel{i0} = char(regexp(odorLabel{i0},sr,'match'));
end

% common odorants used
sharedOdor = {'ethyl acetate','geranyl acetate','ethyl butyrate','isoamyl acetate',...
    '2-heptanone','3-octanol','pentyl acetate','methyl salicylate'};
sharedOr1 = {'Or13a','Or22c','Or24a','Or30a','Or35a','Or42b','Or45a','Or45b',...
    'Or33b-47a','Or49a','Or59a','Or67b','Or74a','Or82a','Or85c'};
sharedOr2 = {'Or13a','Or22c','Or24a','Or30a','Or35a','Or42b','Or45a','Or45b',...
    'Or47a','Or49a','Or59a','Or67b','Or74a','Or82a','Or85c'};

% index of shared odorants in Carlson's data set
inxCarlson = [21,2,22,25,6,18,24,8];
% index of shared odorants in Guangwei's data set
[~,inxGuangwei] = ismember(sharedOdor,odorLabel);%find(contains(odorLabel,sharedOdor));

% index of receptors in Guangwei's data
[~,ixOR1] = ismember(sharedOr1,receptorLabel);   %Guangwei's data
[~,ixOR2] = ismember(sharedOr2,TEXT2);   %Carlson's data

% now we can compare the estimated Sensitivites across different
% odor-receptor paris
sel1 = NUM1(:,ixOR1);
sharedData1 = sel1(inxGuangwei,:);
inxGW = ~isnan(sharedData1);  %index of excitatory interaction

sel2 = NUM2(:,ixOR2);
sharedData2 = sel2(inxCarlson,:);

% corresponding r0
selR0 = NUM2(end,ixOR2);
R0mat = ones(length(sharedOdor),1)*selR0;

% position for strong inhibition
strongInhi = sharedData2 <= -R0mat/2;
% strongInhi = sharedData2 <= -4;

% does the inhibitory interaction overlap with Guangwei's matrix
overLapTest = strongInhi&inxGW;

% estimate the sensitivity matrix from our formula
h = 0.7;
spon = NUM2(end,ixOR2);   %spon taneous activity
spMat = ones(size(sharedData2,1),1)*spon;



Rmax = max(sharedData2(:)) + 1;  %to avoid numerical unstability, so we add 1
alp = Rmax./spon - 1;  % each receptor has one alpha
alpMat = ones(size(sharedData2,1),1)*alp;  %this matrix is used for future use

% ================================================================
% check how the overlap of excitatory interactions in two studies
% ================================================================
thd = 5:5:50;
overlapExci = zeros(length(thd),4);
for i0 = 1:length(thd)
    allInx2 = sharedData2 > thd(i0);
    overlap = allInx2 & inxGW;
    overlapExci(i0,1) = sum(allInx2(:));
    overlapExci(i0,2) = sum(inxGW(:));
    overlapExci(i0,3) = sum(overlap(:))/overlapExci(i0,1);
    overlapExci(i0,4) = sum(overlap(:))/overlapExci(i0,2);
end

% plot how the overlap of excitatory interaction change with parameter
figure
hold on
plot(thd,overlapExci(:,[3,4]),'o-','MarkerSize',12)
lg = legend('eletrophysiology','calcium imaging');
set(lg,'FontSize',16)
xlabel('excitation threshod (Hz)')
ylabel('overlap ratio')
hold off
box on

% number of excitatory interactions in electrophysiology
figure
plot(thd,overlapExci(:,1),'o-','MarkerSize',12)
% lg = legend('eletrophysiology','calcium imaging');
% set(lg,'FontSize',16)
xlabel('excitation threshod (Hz)')
ylabel('counts')
% hold off
box on

% ================================================================
% mismatch of inhibitory interactions in eletrophysiology,but excitaotry in
% calcium imaging method
% ================================================================

% now we need a more reasonable criterion for excitation in
% electrophysiology
exciThd = max(10*ones(size(spMat,1),size(spMat,2)),1.5*spMat);
exciInx = sharedData2 > exciThd;

% now define some probability used to calculate the conditional probability
pEp = sum(exciInx(:))/length(exciInx(:));
pCp = sum(inxGW(:))/length(inxGW(:));
pC0 = 1-pCp;

InhiInx2 = sharedData2 <= -spMat/2; % here is the threshold to define inhibition
pEm = sum(InhiInx2(:))/length(InhiInx2(:));
pE0 = 1-(pEm+pEp);

% now joint probability
temp = inxGW & exciInx;
pEpCp = sum(temp(:))/length(temp(:));
temp = InhiInx2 & inxGW;
pEmCp = sum(temp(:))/length(temp(:));

% now conditional probability
pEm_Cp = pEmCp/pCp;
pEp_Cp = pEpCp/pCp;
pE0_Cp = 1 - pEm_Cp - pEp_Cp;

pEp_C0 = (pEp - pEpCp)/(1-pCp);
pEm_C0 = (pEm - pEmCp)/(1-pCp);
pE0_C0 = 1- pEp_C0 -pEm_C0;

pCp_Ep = pEpCp/pEp;
pC0_Ep = 1- pCp_Ep;

pCp_Em = pEmCp/pEm;
pC0_Em = 1- pCp_Em;

temp = inxGW & ~(InhiInx2|exciInx);
pCpE0 = sum(temp(:))/length(temp(:));
pCp_E0 = pCpE0/pE0;
pC0_E0 = 1 - pCpE0;

% plot the figure

misMatch = InhiInx2 & inxGW;
misRatio = sum(misMatch(:))/sum(InhiInx2(:));

misData = [21,23;7/21,9/23]; % mismatch data
figure
bar(misData(1,:))
set(gca,'XTickLabel',{'1/2','2/3'})
xlabel('threshold ($r_0$)', 'Interpreter','latex')
ylabel('inhibitory interactions')

figure
bar(misData(2,:))
set(gca,'XTickLabel',{'1/2','2/3'})
xlabel('threshold ($r_0$)', 'Interpreter','latex')
ylabel('fraction in elect.')


% reported as zero, but strong inhibition in electrophysiology
InhiInx2 = sharedData2 <= -spMat/3;
misMatch2 = InhiInx2 & ~inxGW;
misRatio = sum(misMatch(:))/sum(InhiInx2(:))
sum(InhiInx2(:))
sum(misMatch2(:))
sum(misMatch(:))/sum(~inxGW(:))

misData = [11/21,12/23]; % mismatch data
figure
bar(misData(1,:))
set(gca,'XTickLabel',{'1/2','2/3'})
xlabel('threshold ($r_0$)', 'Interpreter','latex')
ylabel('faction in elect.')

figure
bar(misData(2,:))
set(gca,'XTickLabel',{'1/2','2/3'})
xlabel('threshold ($r_0$)', 'Interpreter','latex')
ylabel('mismatch ratio')

allInx2 = sharedData2 > 50;
allExci2 = sharedData2(allInx2);
exciW2 = (alpMat(allInx2)./(Rmax./sharedData2(allInx2) - 1) - 1).^(1/h);

InhiInx2 = sharedData2 <= spMat/2;
allInhi2 = max(sharedData2(InhiInx2),0.1);   % for stability
% allInhi = allM(allInx);   % for stability
inhiW2 = ((Rmax./allInhi2-1)./alpMat(InhiInx2) - 1).^(1/h);


% convert into the matrix
NewExciMatr = nan(size(sharedData2,1),size(sharedData2,2));
NewExciMatr(allInx2) = log10(exciW2);

% overlap of the excitatory interactions in two studies
overlapExci = allInx2 & inxGW;
sum(overlapExci(:))/sum(allInx2(:))
sum(overlapExci(:))/sum(inxGW(:))



% Heat map of Wij from Guangwei's data
figure
imagesc(-sharedData1,[2,8])
% xlabel('Or')
set(gca,'XTick',1:size(sharedData1,2),'FontSize',16);
set(gca,'XTickLabel',sharedOr2);
set(gca,'xaxisLocation','top');
set(gca,'YTick',1:size(sharedData1,1),'FontSize',16);
set(gca,'YTickLabel',sharedOdor);
set(gca, 'XTickLabelRotation', 45);
colormap(jet);
c = colorbar; 
c.TickLabels{1} = 'NaN'; 
c.Label.String = 'log10(w)';
c.FontSize = 16;

% Heat map of Wij from Carlson's data
figure
imagesc(NewExciMatr,[0,4])
set(gca,'XTick',1:size(sharedData1,2),'FontSize',16);
set(gca,'XTickLabel',sharedOr2);
set(gca,'xaxisLocation','top');
set(gca,'YTick',1:size(sharedData1,1),'FontSize',16);
set(gca,'YTickLabel',sharedOdor);
set(gca, 'XTickLabelRotation', 45);
colormap(jet);
c = colorbar; 
c.TickLabels{1} = 'NaN'; 
c.Label.String = 'log10(w)';
c.FontSize = 16;

% =========================================================================
% comparing the overlap in these two studies, a scatter plot
% =========================================================================
figure
plot(NewExciMatr(:),-sharedData1(:),'o','MarkerSize',10,'Color',Bu(9,:))
% lg = legend('$$\ln(w)$$')
% set(lg,'Interpreter','latex')
xlabel('electrophysiology')
ylabel('calcium imaging')

% =========================================================================
% Compare spontaneous activity with two diffent measurements
% =========================================================================
% prepare the data
data = [NUM2(end,:)',NUM4(end,:)',sem1(end,:)',sem2(end,:)'];
figure
hold on
errorxy(data,'ColXe',3,'ColYe',4,'MarkSize',10,'EdgeColor',Bu(9,:),'FaceColor',Bu(9,:),...
    'WidthEB',1)
pbaspect([1 1 1])
xL = xlim;
plot(xL',xL','Color',Gr(8,:))
% b1 = NUM2(end,:)'\NUM4(end,:)';
% X = 0:0.2:20;
% plot(X,b1*X,'r')
xlabel('$r_0$ in $10^{-2}$','Interpreter','latex')
ylabel('$r_0$ in $10^{-4}$','Interpreter','latex')


% =========================================================================
% Estimate the sensitivity using data from two concentration
% =========================================================================
thd = 10;  % a simple threshold for the highest concentration
inx1 = NUM2 > thd;
inx2 = NUM4 > 0;  % need to be careful with this threshold
selInx = inx1 & inx2;

% dataVec = [NUM4(selInx),NUM2(selInx)]; %data at two different concentration
rawSel = [NUM4(selInx),NUM2(selInx)];
finalInx = rawSel(:,2) - rawSel(:,1) > 5;
dataVec = rawSel(finalInx,:);  % final selected strong excitation
concentration = ones(length(dataVec),1)*[10^(-4), 10^(-2)];

% here we simply fit the excitatory interaction with hill function with
% hill coeficient of h = 0.7
h = 0.7; 
Rmax = max(dataVec(:))+0.1;
allW = zeros(length(dataVec),1);
X = -10:0.05:-3;
% myFunc = @(x) Rmax*(exp(w)*x).^h./(1+(exp(w)*x).^h);
% f = fittype('R*(exp(w)*x).^0.7./(1+(exp(w)*x).^0.7)');
f = fittype(@(R,w,x) R*(exp(w)*x).^0.7./(1+(exp(w)*x).^0.7));

figure
hold on
numSel = 4;  % only plot 4 lines
plotSelec = randperm(100,4);
allColor = {'red','blue','green','yellow'};
count = 1;   % for color index
for i0 = 1:size(dataVec,1)
    [fit1,gof,fitinfo] = fit(concentration(i0,:)',dataVec(i0,:)',f,'StartPoint',[200,0]);
    allW(i0) = coeffvalues(fit1);
    w = allW(i0);
    myFunc = @(x) Rmax*exp(w+x).^h./(1+(exp(w+x)).^h);
    if  ismember(i0,plotSelec)
        plot(log(concentration(i0,:)'),dataVec(i0,:)','o','MarkerSize',10,...
            'Color',allColor{count},'MarkerFaceColor',allColor{count})
        plot(X,myFunc(X),'LineWidth',2,'Color',allColor{count})
        count = count + 1;
    end
end
hold off
box on
xlabel('$\ln(c)$','Interpreter','latex')
ylabel('firing rate (Hz)')


% histogram of the final excitatory interactions
figure
histogram(allW,15)
legend('h  = 1')
xlabel('$\ln(c)$','Interpreter','latex')
ylabel('Counts')

% ==================================================
% fit both maximum response R and W
% ==================================================
h = 0.7;
f = fittype(@(R,w,x) R*(exp(w)*x).^0.7./(1+(exp(w)*x).^0.7));

figure
hold on
numSel = 4;  % only plot 4 lines
plotSelec = randperm(100,4);
allColor = {'red','blue','green','yellow'};
count = 1;   % for color index
allW = zeros(length(dataVec),2);
for i0 = 1:size(dataVec,1)
    [fit1,gof,fitinfo] = fit(concentration(i0,:)',dataVec(i0,:)',f,...
        'Lower',[50,-10],'Upper',[300,10],'StartPoint',[200,0]);
    allW(i0,:) = coeffvalues(fit1);
    w = allW(i0,2);
    Rmax = allW(i0,1);
    myFunc = @(x) Rmax*exp(w+x).^h./(1+(exp(w+x)).^h);
    if  ismember(i0,plotSelec)
        plot(log(concentration(i0,:)'),dataVec(i0,:)','o','MarkerSize',10,...
            'Color',allColor{count},'MarkerFaceColor',allColor{count})
        plot(X,myFunc(X),'LineWidth',2,'Color',allColor{count})
        count = count + 1;
    end
end
hold off
box on
xlabel('$\ln(c)$','Interpreter','latex')
ylabel('firing rate (Hz)')


% histogram of the final excitatory interactions
figure
histogram(allW(:,2),15)
legend('h  = 0.7')
xlabel('$\ln(c)$','Interpreter','latex')
ylabel('Counts')


% histogram of all Rmax
figure
histogram(allW(:,1),15)
legend('h  = 0.7')
xlabel('$R_{max}$ (Hz)','Interpreter','latex')
ylabel('Counts')


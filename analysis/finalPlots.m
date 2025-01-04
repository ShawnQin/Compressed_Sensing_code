% this program plots all the figures used in the slides. These figures will
% be also used in my PhD thesis

close all
clear

%% define data folder, figure folders and graphics settings
% this section should be ran before the following sections

dFolder = '/Users/shan/Dropbox/olfactionProject/data/gcmi';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% graphics settings
defaultGraphicsSetttings

%define some colors using brewermap
RdBu = brewermap(11,'RdBu');   % red and blue
Bu = brewermap(11,'Blues');    % blues
Gr = brewermap(11,'Greys');    % greys
Rd = brewermap(11,'Reds');
seqColor = brewermap(21,'YlGnBu');

% define my own colors
CN = 255;
myRed = [202,39,44]/CN;
myBrickRed = [242,65,95]/CN;
myOrang = [224,135,51]/CN;
myYellow = [255,241,0]/CN;
myGreen = [10,188,102]/CN;
myCyan = [84,195,241]/CN;
myLightPurple = [232,82,15]/CN;
myDarkPurple = [164,11,93]/CN;
myBlue = [3,110,184]/CN;
myBlueGreen = [76,216,196]/CN;

lBu = [96,166,223]/255; %light blue
dpBu = [63,114,183]/255; % deep blue
dkBu = [50,78,147]/255;   %dark blue
Or = [220,150,71]/255;  % orange
brickRd = [201,69,89]/255;  %brick red
green = [107,169,81]/255;  %green
purple = [113,49,119]/255;  % purple

%% sparisity of optimal W, to show the randomness of W, data from Qianyi's simulation
% this part demonstrates that optimal W is sparse, here we need to select
% some data set, plot the heat map and make a animation of the searching
% process

% used color
% Bu = [16,159,255]/256;
Bu = [143,195,237]/256;
Or = [243,152,0]/256;
% dFolder = '/Users/shan/Dropbox/olfactionProject/data/gcmi/gcmi_N100M20_h1';
dFolder = '/Users/shan/Dropbox/olfactionProject/data/GcmiN100R30-10-03';
% dFolder = '/Users/shan/Dropbox/olfactionProject/data/newGcmi20180831/GcmiSdp0828';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% fName = 'gcmi_whole_N100_R20_S3_sig2_2018-05-22.mat';
fName = 'N100_R30_S2_sig2_2018-10-03.mat';
% fName = 'N100_R20_S3_sig2_2018-08-29.mat';
load(fullfile(dFolder,fName));

% the parameters
N = 100;
M = 30;
sp = 2;
sig =2;

repeats = 100;  %number of permuation

% select one matrix and show its sparsity
w = reshape(allMat(:,1),[M,N]);


% truncates w to be in the range [-3*sigma, 3*sigma]
trucate = true;
if trucate
   trc = 'truc'; %trucation or not
   w(w < -4*sig) = -4*sig;
   w(w > 4*sig) = 4*sig;
else
    trc = 'NOtruc'; %trucation or not
end

% ====================================================
% plot the heatmap of optimal w
% ====================================================
figure
set(gcf,'renderer','Painters')
imagesc(w,[-4,4])
set(gca,'TickLength',[0,0])

% set(gca,'FontSize',16)
% set(gca, 'CLim', [-4, 4]);
% set(gca,'XTick',1:size(w,2),'FontSize',16);
% set(gca,'XTickLabel',[]);
% set(gca,'xaxisLocation','top');
% set(gca,'YTick',1:size(w,1),'FontSize',16);
% set(gca,'YTickLabel',[]);
% set(gca, 'XTickLabelRotation', 45);
colormap(jet);
c = colorbar; 
% c.TickLabels{1} = 'NaN'; 
% c.Label.String = '-log10(EC50)';
% c.FontSize = 16;
xlabel('Odorant')
ylabel('Receptor')
prefix = ['gcmi_exampleW_N',num2str(N),'M',num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ====================================================
% plot the histogram of optimal w
% ====================================================
% for large W, w only need one
% temp = allMat(:);
% hh = histogram(allW,15,'Normalization','pdf',...
%     'EdgeColor','none','FaceColor',lBu,'FaceAlpha',0.4);
% stairs([hh.BinEdges(1),hh.BinEdges,hh.BinEdges(end)],...
%     [0,hh.Values,hh.Values(end),0],'Color',dkBu,'LineWidth',2)
figure
hold on
set(gcf,'renderer','Painters')
h1 = histogram(w(:),60,'Normalization','probability');
h1.FaceColor = lBu; h1.FaceAlpha = 0.4; h1.EdgeColor = 'none';
stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
    [0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',2)

box on
hold off
set(gca,'Layer','top')
xlim([-120,10])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('probability')
prefix = ['gcmi_exampleAllW_N',num2str(N),'M',num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% enlarged
figure
set(gcf,'renderer','Painters')
hold on
h2 = histogram(w(abs(w)<5),20,'Normalization','pdf');
h2.FaceColor = lBu; h2.FaceAlpha = 0.4; h2.EdgeColor = 'none';
stairs([h2.BinEdges(1),h2.BinEdges,h2.BinEdges(end)],...
    [0,h2.Values,h2.Values(end),0],'Color',dpBu,'LineWidth',2)

% get the normal fit parameter
pd = fitdist(w(abs(w)<6),'normal');
X = -6:0.05:6;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y,'Color',Or,'LineWidth',4)
hold off
box on
legend('active w','Gassian fit')
legend boxoff
set(gca,'XLim',[-6,6])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
prefix = ['gcmi_exampleActiveW_N',num2str(N),'M',num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% plot the correlation, row-wise
t1 = corr(w');
y1 = triu(t1,1);
figure
hold on
h3 = histogram(y1(y1~=0),'Normalization','probability');
h3.FaceColor = brickRd; h3.FaceAlpha = 0.4; h3.EdgeColor = 'none';
stairs([h3.BinEdges(1),h3.BinEdges,h3.BinEdges(end)],...
    [0,h3.Values,h3.Values(end),0],'Color',brickRd,'LineWidth',2)
xlim([-1,1])
set(gca,'XTick',-1:0.5:1)
xlabel('correlation among receptors')
ylabel('probability')
prefix = ['gcmi_exampleCorrRecp_N',num2str(N),'M',num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% plot the correlation, colum-wise
t2 = corr(w);
y2 = triu(t2,1);
figure
histogram(y2(y2~=0),'Normalization','probability')
xlim([-1,1])
set(gca,'XTick',-1:0.5:1)
xlabel('correlation among odorants')
ylabel('probability')
prefix = ['gcmi_exampleCorrOdor_N',num2str(N),'M',num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-painters','-depsc',[saveFolder,filesep,prefix,'.eps'])

% merge the two histogram together
figure
histogram(y1(y1~=0),'Normalization','probability','DisplayStyle','stairs','LineWidth',2)
hold on
histogram(y2(y2~=0),'Normalization','probability','DisplayStyle','stairs','LineWidth',2)
hold off
lg = legend('column-wise','row-wise');
set(lg,'FontSize',16)
legend boxoff
xlabel('correlation coefficient')
ylabel('probability')
xlim([-1,1])
set(gca,'XTick',-1:0.5:1)
prefix = ['gcmi_exampleCorrMerge_N',num2str(N),'M',num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-painters','-depsc',[saveFolder,filesep,prefix,'.eps'])


%==================================================================
%correlation coefficients after shuffling, row-wise
%=================================================================

corrType = 'Kendall';  % or Kendall

% first, compare original and shuffuled correlation
L = size(w,1)*size(w,2);

t1 = corr(w','type',corrType);
y1 = triu(t1,1);
rowMean = mean(y1(y1~=0));
rowStd = std(y1(y1~=0));

newW = reshape(w(randperm(L)),size(w,1),size(w,2));
t2 = corr(newW','type',corrType);
y2 = triu(t2,1);

[bc1, edg1] = histcounts(y1(y1~=0),15);
[bc2,edg2] = histcounts(y2(y2~=0),15);

figure
hold on
Y1 = [bc1;bc1]/sum(bc1)/(edg1(2)- edg1(1)); %normalized, probability
ar1 = area(sort(edg1([1:end-1 2:end])),Y1(1:end));
ar1.FaceAlpha = 0.3;
ar1.FaceColor = lBu;
ar1.EdgeColor = dpBu;
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
xlabel('$\tau$','Interpreter','latex')
ylabel('pdf')
% xlim([-1,1])
set(gca,'XTick',-1:0.5:1)
% h5 = histogram(y2(y2~=0),'Normalization','probability','DisplayStyle','bar',...
%     'EdgeColor','none','FaceColor',Gr(7,:));
% h5.FaceAlpha = 0.3;
% stairs([h5.BinEdges(1),h5.BinEdges,h5.BinEdges(end)],...
%     [0,h5.Values,h5.Values(end),0],'Color',Gr(7,:),'LineWidth',2)
% 
% h4 = histogram(y1(y1~=0),'Normalization','probability','DisplayStyle','bar',...
%     'EdgeColor','none');
% h4.FaceAlpha = 0.4; h4.FaceColor = lBu;
% stairs([h4.BinEdges(1),h4.BinEdges,h4.BinEdges(end)],...
%     [0,h4.Values,h4.Values(end),0],'Color',dpBu,'LineWidth',2)
% xlim([-1,1])
% 
% hold off
% lg = legend('shuffled','original');
% set(lg,'FontSize',16)
% legend boxoff
% box on
% xlabel('$\rho_s$','Interpreter','latex')
% ylabel('probability')
% xlim([-1,1])
% set(gca,'XTick',-1:0.5:1)
prefix = ['gcmi_exampleCorrComp_RowShuff_',corrType,'_',trc,'_N',...
    num2str(N),'M',num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-painters','-dpdf',[saveFolder,filesep,prefix,'.pdf'])


% statistics after shuffling
repeats = 100;
% L = size(w,1)*size(w,2);
allMeanCorr = zeros(repeats,1);
allStdCorr = zeros(repeats,1);

for i0 = 1:repeats
    newW = reshape(w(randperm(L)),size(w,1),size(w,2));
    temp = corr(newW','type',corrType);
    y3 = triu(temp,1);
    allMeanCorr(i0) = mean(y3(y3~=0));
    allStdCorr(i0) = std(y3(y3~=0));
end
figure
hold on
scatter(allMeanCorr,allStdCorr,20,[0.5,0.5,0.5],'filled')
xlim([-0.05,0.05])
ylim([0.04,0.11])
set(gca,'XTick',-0.05:0.05:0.05)
scatter(rowMean,rowStd,80,dkBu,'filled')
hold off
box on
lg = legend('shuffled','original');
xlabel('$\langle \rho_s \rangle$','Interpreter','latex')
ylabel('$\sigma_{\rho_s}$','Interpreter','latex')
prefix = ['gcmi_exampleCorrComp_RowShuffScatter_',corrType,'_',trc,'_N',...
    num2str(N),'M',num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-painters','-depsc',[saveFolder,filesep,prefix,'.eps'])

%==================================================================
%correlation coefficients after shuffling, cloumn-wise
%=================================================================
% first, compare original and shuffuled correlation
% corrType = 'Spearman';  % or Spearman
t1 = corr(w,'type',corrType);
y1 = triu(t1,1);
rowMean = mean(y1(y1~=0));
rowStd = std(y1(y1~=0));

newW = reshape(w(randperm(L)),size(w,1),size(w,2));
t2 = corr(newW,'type',corrType);
y2 = triu(t2,1);

[bc1, edg1] = histcounts(y1(y1~=0),25);
[bc2,edg2] = histcounts(y2(y2~=0),25);

figure
hold on

Y1 = [bc1;bc1]/sum(bc1)/(edg1(2)- edg1(1)); %normalized, probability
ar1 = area(sort(edg1([1:end-1 2:end])),Y1(1:end));
ar1.FaceAlpha = 0.3;
ar1.FaceColor = lBu;
ar1.EdgeColor = dpBu;
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
xlabel('$\tau$','Interpreter','latex')
ylabel('pdf')
% xlim([-1,1])
set(gca,'XTick',-1:0.5:1)
hold off
% h2 = histogram(y2(y2~=0),'Normalization','probability','DisplayStyle','bar',...
%     'EdgeColor','none','FaceColor',Gr(7,:),'FaceAlpha',0.3);
% stairs([h2.BinEdges(1),h2.BinEdges,h2.BinEdges(end)],...
%     [0,h2.Values,h2.Values(end),0],'Color',Gr(7,:),'LineWidth',2)
% 
% h1 = histogram(y1(y1~=0),'Normalization','probability','DisplayStyle','bar',...
%     'EdgeColor','none','FaceColor',lBu,'FaceAlpha',0.4);
% stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
%     [0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',2)
% xlim([-1,1])
% 
% hold off
% box on
% lg = legend('shuffled','original');
% set(lg,'FontSize',16)
% legend boxoff
% xlabel('$\rho_s$','Interpreter','latex')
% ylabel('probability')
% xlim([-1,1])
% set(gca,'XTick',-1:0.5:1)
prefix = ['gcmi_CorrComp_columnShuff_',corrType,'_',trc,'_N',num2str(N),'M',...
    num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-painters','-depsc',[saveFolder,filesep,prefix,'.eps'])
print('-dpdf',[saveFolder,filesep,prefix,'.pdf'])


% statistics after shuffling
repeats = 100;
L = size(w,1)*size(w,2);
allMeanCorr = zeros(repeats,1);
allStdCorr = zeros(repeats,1);

for i0 = 1:repeats
    newW = reshape(w(randperm(L)),size(w,1),size(w,2));
    temp = corr(newW,'type',corrType);
    y3 = triu(temp,1);
    allMeanCorr(i0) = mean(y3(y3~=0));
    allStdCorr(i0) = std(y3(y3~=0));
end
figure
hold on
box on
scatter(allMeanCorr,allStdCorr,20,[0.5,0.5,0.5],'filled')
xlim([-0.01,0.01])
ylim([0.12,0.18])
set(gca,'XTick',-0.01:0.01:0.01)
scatter(rowMean,rowStd,80,dkBu,'filled')
hold off
lg = legend('shuffled','original');
xlabel('$\langle \rho_s \rangle$','Interpreter','latex')
ylabel('$\sigma_{\rho_s}$','Interpreter','latex')
prefix = ['gcmi_exampleCorrComp_ColumnShuffScatter_',corrType,'_',trc,'_N',...
    num2str(N),'M',num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% ========================================================
% Compare the eigenvalues of optimal and permutated matrix
% Row-wise
% ========================================================
repeats = 500;  %number of permuation

% select one matrix and show its sparsity
w = reshape(allMat(:,1),[M,N]);

% truncates w to be in the range [-4*sigma, 4*sigma]
trucate = true;   %modified on 12/06/2018
if trucate
   trc = 'truc'; %trucation or not
   w(w < -4*sig) = -4*sig;
   w(w > 4*sig) = 4*sig;
else
    trc = 'NOtruc'; %trucation or not
end

corrType = 'Kendall';  %or Spearman

C1 = corr(w','type',corrType);
eigenOpt = eig(C1);

% random permutation
eigenPerm = zeros(M, repeats);
aveOpt = [mean(eigenOpt),std(eigenOpt)];

for i0 = 1:repeats
    newW = reshape(w(randperm(M*N)),size(w,1),size(w,2));
    temp = corr(newW','type',corrType);
    eigenPerm(:,i0) = eig(temp);
end
avePerm = [mean(eigenPerm,1);std(eigenPerm,0,1)];

figure
set(gcf,'Renderer','painters')
hold on
histogram(avePerm(2,:),'FaceColor',Bu,'Normalization','probability','LineWidth',0.5);
ah = gca;
ylim = ah.YLim;
plot([aveOpt(2); aveOpt(2)],ylim,'r--','LineWidth',2)
box on
set(gca,'XLim',[0.25, 0.5])
legend('shuffle','original')
xlabel('$\sigma_{\lambda}$','Interpreter','latex')
ylabel('probability')
prefix = ['gcmi_CorrComp_Row_eigen_',corrType,'_',trc,'_N',num2str(N),'M',...
    num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% CDF compare
figure
hold on
for i0 = 1:repeats
    [f,x] = ecdf(eigenPerm(:,i0));
    plot(x,f,'LineWidth',1,'Color',Gr(6,:))
end

[f,x] = ecdf(eigenOpt);
plot(x,f,'LineWidth',2,'Color','red')
hold off
lg = legend('shuffle','original');
legend boxoff
xlabel('$\lambda $','Interpreter','latex')
ylabel('ECDF')
prefix = ['gcmi_CorrComp_Row_CDF_',corrType,'_',trc,'_N',num2str(N),'M',...
    num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ==================================================
% column-wise correlation matrix, eigenvalues comparison
% ==================================================
thd = 1e-5;  % threshold used to determine if a eigenvalue is zero
C1 = corr(w,'type',corrType);
eigenOpt = eig(C1);

% random permutation
eigenPerm = zeros(N, repeats);
aveOpt = [mean(eigenOpt(72:end)),std(eigenOpt(72:end))];

for i0 = 1:repeats
    newW = reshape(w(randperm(M*N)),size(w,1),size(w,2));
    temp = corr(newW,'type',corrType);
    eigenPerm(:,i0) = eig(temp);
end
avePerm = [mean(eigenPerm,1);std(eigenPerm,0,1)];

figure
set(gcf,'Renderer','painters')
hold on
histogram(avePerm(2,72:end),'FaceColor',Bu,'Normalization','probability','LineWidth',0.5);
ah = gca;
ylim = ah.YLim;
plot([aveOpt(2); aveOpt(2)],ylim,'r--','LineWidth',2)
box on
set(gca,'XLim',[0.9, 1.6])
legend('shuffle','original')
xlabel('$\sigma_{\lambda}$','Interpreter','latex')
ylabel('probability')
prefix = ['gcmi_CorrComp_Column_eigen_',corrType,'_',trc,'_N',num2str(N),...
    'M',num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


figure
set(gcf,'Renderer','painters')
hold on
for i0 = 1:repeats
    [f,x] = ecdf(eigenPerm(72:end,i0));
    plot(x,f,'LineWidth',1,'Color',Gr(6,:))
end

[f,x] = ecdf(eigenOpt(72:end));
plot(x,f,'LineWidth',2,'Color','red')
hold off
lg = legend('shuffle','original');
legend boxoff
xlabel('$\lambda $','Interpreter','latex')
ylabel('ECDF')
prefix = ['gcmi_CorrComp_Column_CDF_',corrType,'_',trc,'_N',num2str(N),'M',...
    num2str(M),'sig',num2str(sig),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])



% ==============================================
% map the elements of original matrix into a stardard 
% Gaussian distribution and compare the eigenvalues
% ==============================================
[f,x] = ecdf(w(:));
f(end) = 1-1e-4;
[~, inx] = sort(w(:));   %index
newW = norminv(f(2:end));
w0 = zeros(size(w,1),size(w,2));
w0(inx) = newW;

C0 = corr(w0','type',corrType);
eigenOpt = eig(C0);
eigenPerm = zeros(M, repeats);
aveOpt = [mean(eigenOpt),std(eigenOpt)];

for i0 = 1:repeats
    newW = reshape(w0(randperm(M*N)),size(w0,1),size(w0,2));
    temp = corr(newW','type',corrType);
    eigenPerm(:,i0) = eig(temp);
end
avePerm = [mean(eigenPerm,1);std(eigenPerm,0,1)];
largest = max(eigenPerm);  %largest eigen values

% compare the standard deviation of eigenvalues
figure
hold on
histogram(avePerm(2,:),'FaceColor',Bu,'Normalization','probability');
ah = gca;
ylim = ah.YLim;
plot([aveOpt(2); aveOpt(2)],ylim,'r--','LineWidth',2)
box on
legend('shuffle','original')
xlabel('$\sigma_{\lambda}$','Interpreter','latex')
ylabel('probability')

% compare the maxium eigenvalue
figure
hold on
histogram(largest,'FaceColor',Bu,'Normalization','probability');
ah = gca;
ylim = ah.YLim;
plot([max(eigenOpt); max(eigenOpt)],ylim,'r--','LineWidth',2)
box on
legend('shuffle','original')
xlabel('$\lambda_{max}$','Interpreter','latex')
ylabel('probability')


% ecdf of data and Gaussian
figure
plot(x,f)
xlabel('w')
ylabel('ecdf')

% standard Gaussian
pd = makedist('Normal');
X = -3:.1:3;
cdf_normal = cdf(pd,X);
figure
plot(X,cdf_normal)
xlabel('x')
ylabel('cdf')

%% How sparsity of W change with N
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% load the summary data set
dName = 'gcmi_Ndp_M20Sp2sig2_12-Sep-2018.mat';
load(fullfile(dFolder,dName))

M = 20;
sp = 2;
sig = 2;


labelFontSize = 24;

odorInx = 1:1:8; %odor selected
colorInx = 9;
LineWidth = 2;
% sp of W vs sigmac
figure
errorbar(allN(odorInx)',meanSpW(odorInx),stdSpW(odorInx),'o-','MarkerSize',12,...
        'MarkerFaceColor',Bu(colorInx,:),'Color',Bu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend(['M=',num2str(M),',n=',num2str(3),',\sigma_c =',num2str(2)],'Location','southwest');
legend boxoff
xlabel('$N$','interpreter','latex')
ylabel('$\rho_{w}$','Interpreter','latex')
prefix = ['gcmiNdp_M',num2str(M),'sp',num2str(sp),'sig',num2str(sig),'_rhoW',date];
saveas(gcf,[sFolder,filesep,prefix,'.fig'])
print('-depsc',[sFolder,filesep,prefix,'.eps'])

% set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

% mean W vs sigmac
figure
errorbar(allN(odorInx)',meanW(odorInx),stdW(odorInx),'o-','MarkerSize',12,...
        'MarkerFaceColor',Bu(colorInx,:),'Color',Bu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend(['M=',num2str(M),',n=',num2str(3),',\sigma_c =',num2str(2)],'Location','southwest');
legend boxoff
xlabel('$N$','interpreter','latex')
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
prefix = ['gcmiNdp_M',num2str(M),'sp',num2str(sp),'sig',num2str(sig),'_muW',date];
saveas(gcf,[sFolder,filesep,prefix,'.fig'])
print('-depsc',[sFolder,filesep,prefix,'.eps'])

% set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

%std w vs sigmac
figure
errorbar(allN(odorInx)',meanSigW(odorInx),stdSigW(odorInx),'o-','MarkerSize',12,...
        'MarkerFaceColor',Bu(colorInx,:),'Color',Bu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend(['M=',num2str(M),',n=',num2str(3),',\sigma_c =',num2str(2)],'Location','southwest');
legend boxoff
xlabel('$N$','Interpreter','latex')
ylabel('$\sigma_w$','Interpreter','latex')
% ylabel('mean W','FontSize',labelFontSize)
% set(gca,'FontSize',20,'LineWidth',1.5,'XScale','linear')

%save the figure
prefix = ['gcmiNdp_M',num2str(M),'sp',num2str(sp),'sig',num2str(sig),'_sigW',date];
saveas(gcf,[sFolder,filesep,prefix,'.fig'])
print('-depsc',[sFolder,filesep,prefix,'.eps'])

%% How differential entropy change with N, sigma_c, using Qianyi's data
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% load the summary data set
% dName = 'gcmi_diffN_diffSig_M30Summary_QY_10-Aug-2018.mat';
dName = 'N50M13sp2_Summary_diffSig_2018-10-05.mat';
load(fullfile(dFolder,dName))

M = 13;
sp = 2;
sig = 2;

% ================================
% changes with N
% ================================
ix = 1;  % select sigma = 2
figure
errorbar(allN,allMeanFmin(4,:),allStdFmin(4,:),'o-','MarkerSize',12,'Color',Bu(10,:), ...
    'MarkerFaceColor',Bu(9,:),'LineWidth',2,'CapSize',0)
lg = legend('M=30,sp=2,\sigma_c = 2','Location','southeast');
set(lg,'FontSize',20)
legend boxoff
xlabel('$N$','interpreter','latex')
ylabel('$f_{min}$','interpreter','latex')
prefix = ['Gcmi_fmin_Ndp_M',num2str(M),'sig',num2str(sig),'_',date,'.mat'];
print(gcf,'-depsc',fullfile(sFolder,[prefix,'.eps']))
saveas(gcf,fullfile(sFolder,[prefix,'.fig']))


% ================================
% changes with sig
% ================================
ixN = 5;  % select sigma = 2
figure
errorbar(allSig',allMeanFmin(:,ixN),allStdFmin(:,ixN),'o-','MarkerSize',12,'Color',Bu(10,:), ...
    'MarkerFaceColor',Bu(9,:),'LineWidth',2,'CapSize',0)
lg = legend('M=30,sp=2,\sigma_c = 2','Location','southeast');
set(lg,'FontSize',20)
legend boxoff
xlabel('$\sigma_c$','interpreter','latex')
ylabel('$f_{min}$','interpreter','latex')
prefix = ['Gcmi_fmin_Sigdp_M',num2str(M),'sig',num2str(sig),'_',date,'.mat'];
print(gcf,'-depsc',fullfile(sFolder,[prefix,'.eps']))
saveas(gcf,fullfile(sFolder,[prefix,'.fig']))

%% How sparsity of W changes with input sigma_c, using simulation data, the gcmi algorithm
% load the summarized data
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
fName = 'gcmi_diffN_summary_06222018.mat';

load(fullfile(dFolder,fName));

% only select part of the data to show
inx = 4:6;
% ========================================================================
% plot how sparsity of optimal W changes with sigma_c with selected N
% ========================================================================
figure
hold on
for i0 = 1:length(inx)
    errorbar(allSig,meanSp(:,inx(i0)),stdSp(:,inx(i0)),'o-','MarkerSize',...
        12,'MarkerFaceColor',Bu(2+3*i0,:),'Color',Bu(2+3*i0,:),'LineWidth',2,...
    'CapSize',0)
end
hold off
lg = legend(['N=',num2str(nOdor(inx(1)))],['N=',num2str(nOdor(inx(2)))],...
    ['N=',num2str(nOdor(inx(3)))],'Location','northwest');
legend boxoff
box on
set(lg,'FontSize',20)
xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('sparsity of W','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_diffN_sp_sigc_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ========================================================================
% plot how the peak changes with sigmac
% ========================================================================
figure
hold on
for i0 = 1:length(inx)
    errorbar(allSig,meanPeak(:,inx(i0)),stdPeak(:,inx(i0)),'o-','MarkerSize',...
        12,'MarkerFaceColor',Bu(2+3*i0,:),'Color',Bu(2+3*i0,:),'LineWidth',2,...
    'CapSize',0)
end
hold off
lg = legend(['N=',num2str(nOdor(inx(1)))],['N=',num2str(nOdor(inx(2)))],...
    ['N=',num2str(nOdor(inx(3)))],'Location','southwest');
legend boxoff
box on
set(lg,'FontSize',20)
xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('peak of $\ln{w}$','Interpreter','latex','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_diffN_muW_sigc_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ========================================================================
% plot how the standard deviation of W changes with sigmac
% ========================================================================
figure
hold on
for i0 = 1:length(inx)
    errorbar(allSig,meanSig(:,inx(i0)),stdSig(:,inx(i0)),'o-','MarkerSize',...
        12,'MarkerFaceColor',Bu(2+3*i0,:),'Color',Bu(2+3*i0,:),'LineWidth',2,...
        'CapSize',0)
end
hold off
lg = legend(['N=',num2str(nOdor(inx(1)))],['N=',num2str(nOdor(inx(2)))],...
    ['N=',num2str(nOdor(inx(3)))],'Location','northwest');
legend boxoff
box on
set(lg,'FontSize',20)
xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('$\sigma_{\ln{(w)}}$','Interpreter','latex','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_diffN_sigW_sigc_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% now we plot the statistics as a function of N for a give sigmac
inxSig  = 2:2:8;
inxOdor = 4:7;

% ========================================================================
% first plot how the sparsity change with N, for different sigmac
% ========================================================================
figure
hold on
for i0 = 1:length(inxSig)
    errorbar(nOdor(inxOdor),meanPeak(inxSig(i0),inxOdor)',stdPeak(inxSig(i0),inxOdor)','o-','MarkerSize',...
        12,'MarkerFaceColor',Bu(3+inxSig(i0),:),'Color',Bu(3+inxSig(i0),:),...
        'CapSize',0,'LineWidth',2)
end
hold off
lg = legend(['\sigma_c=',num2str((allSig(inxSig(1))))],['\sigma_c=',num2str((allSig(inxSig(2))))],...
    ['\sigma_c=',num2str((allSig(inxSig(3))))],['\sigma_c=',num2str((allSig(inxSig(4))))],...
    'Location','southeast');
legend boxoff
box on
set(lg,'FontSize',16)
xlabel('number of odorants','FontSize',28)
ylabel('peak of active W','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_diffSigc_muW_N_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ========================================================================
% plot how the position of active w changes with N for different sigmac
% ========================================================================
figure
hold on
for i0 = 1:length(inxSig)
    errorbar(nOdor(inxOdor),meanSp(inxSig(i0),inxOdor)',stdSp(inxSig(i0),inxOdor)','o-','MarkerSize',...
        12,'MarkerFaceColor',Bu(3+inxSig(i0),:),'Color',Bu(3+inxSig(i0),:),...
        'CapSize',0,'LineWidth',2)
end
hold off
lg = legend(['\sigma_c=',num2str((allSig(inxSig(1))))],['\sigma_c=',num2str((allSig(inxSig(2))))],...
    ['\sigma_c=',num2str((allSig(inxSig(3))))],['\sigma_c=',num2str((allSig(inxSig(4))))],...
    'Location','southwest');
legend boxoff
box on
set(lg,'FontSize',16)
xlabel('number of odorants','FontSize',28)
ylabel('sparsity of W','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_diffSigc_spW_N_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ========================================================================
% plot how the active elements change with sigmac, for different sigmac
% ========================================================================
figure
hold on
for i0 = 1:length(inxSig)
    errorbar(nOdor(inxOdor),meanSen(inxSig(i0),inxOdor)',stdSen(inxSig(i0),inxOdor)','o-','MarkerSize',...
        12,'MarkerFaceColor',Bu(3+inxSig(i0),:),'Color',Bu(3+inxSig(i0),:),...
        'CapSize',0,'LineWidth',2)
end
hold off
lg = legend(['\sigma_c=',num2str((allSig(inxSig(1))))],['\sigma_c=',num2str((allSig(inxSig(2))))],...
    ['\sigma_c=',num2str((allSig(inxSig(3))))],['\sigma_c=',num2str((allSig(inxSig(4))))],...
    'Location','northwest');
legend boxoff
box on
set(lg,'FontSize',16)
xlabel('number of odorants','FontSize',28)
ylabel({'active receptors', 'per odor'},'FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_diffSigc_activeWperOdor_N_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% ========================================================================
% plot how the standard divation of W change with sigmac, for different sigmac
% ========================================================================
figure
hold on
for i0 = 1:length(inxSig)
    errorbar(nOdor(inxOdor),meanSig(inxSig(i0),inxOdor)',stdSig(inxSig(i0),inxOdor)','o-','MarkerSize',...
        12,'MarkerFaceColor',Bu(3+inxSig(i0),:),'Color',Bu(3+inxSig(i0),:),....
        'CapSize',0,'LineWidth',2)
end
hold off
lg = legend(['\sigma_c=',num2str((allSig(inxSig(1))))],['\sigma_c=',num2str((allSig(inxSig(2))))],...
    ['\sigma_c=',num2str((allSig(inxSig(3))))],['\sigma_c=',num2str((allSig(inxSig(4))))],...
    'Location','northwest');
legend boxoff
ylim([1,2])
box on
set(lg,'FontSize',16)
xlabel('number of odorants','FontSize',28)
ylabel('$\sigma_{\ln{(w)}}$','Interpreter','latex','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_diffSigc_stdW_N_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ========================================================================
% plot how the sparsity of W change with N, for different sigmac, using all
% the data
% ========================================================================
inxSig  = [2,4,7];
inxOdor = 1:7;
figure
hold on
for i0 = 1:length(inxSig)
    errorbar(nOdor(inxOdor),meanSp(inxSig(i0),inxOdor)',stdSp(inxSig(i0),inxOdor)','o-','MarkerSize',...
        12,'MarkerFaceColor',Bu(3+inxSig(i0),:),'Color',Bu(3+inxSig(i0),:),...
        'CapSize',0,'LineWidth',2)
end
hold off
lg = legend(['\sigma_c=',num2str((allSig(inxSig(1))))],['\sigma_c=',num2str((allSig(inxSig(2))))],...
    ['\sigma_c=',num2str((allSig(inxSig(3))))],...
    'Location','southeast');
legend boxoff
box on
set(lg,'FontSize',16)
xlabel('number of odorants','FontSize',28)
ylabel('sparsity of active W','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_diffSigc_spW_N_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ========================================================================
% plot how the std of W change with N, for different sigmac, using all
% the data
% ========================================================================
inxSig  = [2,4,7];
inxOdor = 1:7;
figure
hold on
for i0 = 1:length(inxSig)
    errorbar(nOdor(inxOdor),meanSig(inxSig(i0),inxOdor)',stdSig(inxSig(i0),inxOdor)','o-','MarkerSize',...
        12,'MarkerFaceColor',Bu(3+inxSig(i0),:),'Color',Bu(3+inxSig(i0),:),...
        'CapSize',0,'LineWidth',2)
end
hold off
lg = legend(['\sigma_c=',num2str((allSig(inxSig(1))))],['\sigma_c=',num2str((allSig(inxSig(2))))],...
    ['\sigma_c=',num2str((allSig(inxSig(3))))],...
    'Location','southeast');
legend boxoff
box on
set(lg,'FontSize',16)
xlabel('number of odorants','FontSize',28)
ylabel('sparsity of active W','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_diffSigc_sigW_N_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% How sparisty of W change with input sparsity (n), GCMI algorithm, N=100,M=20,sigma_c = 2

%=============================================================
%prepare the data
%=============================================================
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

fileName = 'N100M20Summary_h1_28-Jun-2018.mat';
load(fullfile(dFolder,fileName))

N = 100;
M = 20;
sig = 2;
h = 1;
% ================================================
% plot how sparisty of W depends on input sparisty
% ================================================
figure
errorbar(dataSumm.sp',dataSumm.spMean,dataSumm.spStd,'o-','MarkerSize',...
      12,'MarkerFaceColor',Bu(10,:),'Color',Bu(10,:),'CapSize',0,'LineWidth',2)
lg = legend(['N=',num2str(N),',M=',num2str(M),',\sigma_c = ',num2str(sig),',h=',num2str(h)]);
legend boxoff
set(lg,'FontSize',20)
xlim([1,10])
box on
xlabel('sparisity of input (n)','FontSize',28)
ylabel('sparisty of W','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
prefix = ['gcmi_spDepen_N',num2str(N),'M',num2str(M),'sig',num2str(sig),'h',...
    num2str(h),'_spW_',date];
saveas(gcf,[outFolder,filesep,prefix,'.fig'])
print('-depsc',[outFolder,filesep,prefix,'.eps'])


% ================================================
% plot how average of active W depends on input sparisty
% ================================================
figure
errorbar(dataSumm.sp',dataSumm.meanAveW,dataSumm.stdAveW,'o-','MarkerSize',...
      12,'MarkerFaceColor',Bu(10,:),'Color',Bu(10,:),'CapSize',0,'LineWidth',2)
lg = legend(['N=',num2str(N),',M=',num2str(M),',\sigma_c = ',num2str(sig),',h=',num2str(h)]);
legend boxoff
set(lg,'FontSize',20)
xlim([1,10])
box on
xlabel('sparisity of input (n)','FontSize',28)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
prefix = ['gcmi_spDepen_N',num2str(N),'M',num2str(M),'sig',num2str(sig),'h',...
    num2str(h),'_muW_',date];
saveas(gcf,[outFolder,filesep,prefix,'.fig'])
print('-depsc',[outFolder,filesep,prefix,'.eps'])


% ================================================
% plot how std of active W depends on input sparisty
% ================================================
figure
errorbar(dataSumm.sp',dataSumm.meanStdW,dataSumm.stdStdW,'o-','MarkerSize',...
      12,'MarkerFaceColor',Bu(10,:),'Color',Bu(10,:),'CapSize',0,'LineWidth',2)
lg = legend(['N=',num2str(N),',M=',num2str(M),',\sigma_c = ',num2str(sig),',h=',num2str(h)]);
legend boxoff
set(lg,'FontSize',20)
xlim([1,10])
box on
xlabel('sparisity of input (n)','FontSize',28)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
prefix = ['gcmi_spDepen_N',num2str(N),'M',num2str(M),'sig',num2str(sig),'h',...
    num2str(h),'_sigW_',date];
saveas(gcf,[outFolder,filesep,prefix,'.fig'])
print('-depsc',[outFolder,filesep,prefix,'.eps'])

%% GCMI with paramterized distribution, search
% this part plot the summary data using simulation with paramerized
% distribution of W
%
close all
% clear

dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';


% =================================================================
% How all the parameters of W changes with the input sigma_c
% =================================================================
fileName = 'gcmi_distri_summData_28-Aug-2018.mat';
load(fullfile(dFolder,fileName))

N = 100;
M = 20;
sp = 2;

% tight plot, 2x2 figures
figure
figureSize = [0 0 11 8];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the marginsre
[ha, pos] = tight_subplot(2,2,[.05 .12],[.1 .05],[.12 .05]);


% sparisty of W
axes(ha(1))
errorbar(summData.allSig',summData.meanRho,summData.stdRho,'o-','MarkerSize',12,...
        'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,...
        'CapSize',0)
lg = legend('N = 100,M= 20, n = 2','Location','north');
set(lg,'FontSize',16)
legend boxoff
% xlabel('$\sigma_c$','interpreter','latex')
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',24)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear','XTick',[])


% mean value of W
axes(ha(2))
errorbar(summData.allSig',summData.meanW,summData.stdW,'o-','MarkerSize',12,...
        'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,...
        'CapSize',0)
% xlabel('$\sigma_c$','interpreter','latex')
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',24)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear','XTick',[])

% std of W
axes(ha(3))
errorbar(summData.allSig',summData.meanSigW,summData.stdSigW,'o-','MarkerSize',12,...
        'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,...
        'CapSize',0)
xlabel('$\sigma_c$','interpreter','latex')
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',24)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')



axes(ha(4))
errorbar(summData.allSig',summData.meanFmin,summData.stdFmin,'o-','MarkerSize',12,...
        'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,...
        'CapSize',0)
xlabel('$\sigma_c$','interpreter','latex')
ylabel('$f_{min}$','Interpreter','latex','FontSize',24)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')


prefix = ['gcmi_distr_sigDpd_N',num2str(N),'M',num2str(M),'sp',num2str(sp),'_',date];
saveas(gcf,[sFolder,filesep,prefix,'.fig'])
print('-depsc',[sFolder,filesep,prefix,'.eps'])

% =================================================================
% How all the parameters of W changes with the input N
% =================================================================

fileName = 'gmci_distri_N_summData_10-Aug-2018.mat';
load(fullfile(dFolder,fileName))

% N = 100;
M = 30;
sp = 5;
sig = 2;

% tight plot, 2x2 figures
figure
figureSize = [0 0 11 8];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(2,2,[.05 .12],[.1 .05],[.12 .05]);


% sparisty of W
axes(ha(1))
errorbar(allN',allMeanSpW,allStdSpW,'o-','MarkerSize',12,...
        'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,...
        'CapSize',0)
lg = legend('M= 30, n = 5,\sigma_c = 2','Location','southeast');
set(lg,'FontSize',16)
legend boxoff
% xlabel('$\sigma_c$','interpreter','latex')
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',24)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear','XTick',[])


% mean value of W
axes(ha(2))
errorbar(allN',allMeanAveW,allStdAveW,'o-','MarkerSize',12,...
        'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,...
        'CapSize',0)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',24)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear','XTick',[])

% std of W
axes(ha(3))
errorbar(allN',allMeanSigW,allStdSigW,'o-','MarkerSize',12,...
        'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,...
        'CapSize',0)
xlabel('$N$','interpreter','latex')
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',24)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')



axes(ha(4))
errorbar(allN',allMeanFmin,allStdFmin,'o-','MarkerSize',12,...
        'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,...
        'CapSize',0)
xlabel('$N$','interpreter','latex')
ylabel('$f_{min}$','Interpreter','latex','FontSize',24)
set(gca,'FontSize',24,'LineWidth',1.5)


prefix = ['gcmi_distr_NDpd_',num2str(N),'M',num2str(M),'sp',num2str(sp),'sig',num2str(sig),'_',date];
saveas(gcf,[sFolder,filesep,prefix,'.fig'])
print('-depsc',[sFolder,filesep,prefix,'.eps'])

%% prepare and save the data
dataFolder = '/Users/shan/Dropbox/olfactionProject/data/gcmi/gcmi_N100M20_h1';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
filesRaw = {allFile.name}';

% screen the files, only part are used. Since the folder may contain
% heterogeneous data set
str_marker = '2018-05-22';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

% define the basic parameters
% R = [60,100,150,200];
R = 20;      % number of receptors
num = 40;    % number of repeat
N = 100;     % number of odorants
sig = 2;     % sigma_c
h = 1;       % Hill coefficient of input function

% data structure store all the data
dataSumm = struct('allFmean',[],'allFstd',[],'sp',[],'allSig',[],'spMean',[],'spStd',[]);
dataSumm.allSig = 2;  % sigma_c
% string flag to extract information
spStr = '(?<= *_S)[\d.]+(?=_)';

% cycle through all the data file, get R and allSig information
allSpInput = [];  %sparsity of input
for i0 = 1:length(files)
    allSpInput = [allSpInput, str2num(char(regexp(files{i0},spStr,'match')))];
end
% [val1,~,Rinx] = unique(allR);
[val2,~,Sinx] = unique(allSpInput);
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
        thd = -5*sig;
        
        sp = sum(temp.allMat > thd,1)/R/N;
%         dataSumm.sparsity(Rinx(i0),Sinx(i0)) = sum(tm > thd)/length(tm);
        dataSumm.spMean(Sinx(i0),1) = mean(sp);
        dataSumm.spStd(Sinx(i0),1) = std(sp);
        
        % save the matrix columwize
        dataSumm.allW{Sinx(i0),1} = temp.allMat;  %normalize

    end
 
end
% save the data
dataName = fullfile(saveFolder,['N100M20Summary_h',num2str(h),'_',date,'.mat']);
% dataName = fullfile(saveFolder,'N100M20Summary_h2_05212018.mat');
save(dataName,'dataSumm')


%================================================
% plot how sparisty of W change with n
%================================================
figure
errorbar(dataSumm.sp',dataSumm.spMean,dataSumm.spStd,'o-','MarkerSize',12,...
    'Color',Bu(10,:),'MarkerFaceColor',Bu(10,:),'LineWidth',2,'CapSize',0)
legend('N=100,M=20,\sigma_c = 2, h = 1')
legend boxoff
xlabel('number of odors (sparsity)')
ylabel('sparsity of W')
figPref = ['N100M20_sparistyW_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])


%================================================
% mean and std of active elements
%================================================
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

% plot mean value of active W
figure
errorbar(dataSumm.sp',allMeanWmean,allMeanWstd,'o-','MarkerSize',12,...
    'Color',Bu(10,:),'MarkerFaceColor',Bu(10,:),'LineWidth',2,'CapSize',0)
legend('N=100,M=20,\sigma_c = 2, h = 1')
legend boxoff
xlabel('number of odors (sparsity)')
ylabel('$\mu_w$','Interpreter','latex')
figPref = ['N100M20_meanActiveW_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

% plot the std of W
figure
errorbar(dataSumm.sp',allStdWmean,allStdWstd,'o-','MarkerSize',12,'LineWidth',2,...
    'Color',Bu(10,:),'MarkerFaceColor',Bu(10,:),'CapSize',0)
ylim([1.5,3])
legend('N=100,M=20,\sigma_c = 2, h = 1')
legend boxoff
xlabel('number of odors (sparsity)')
ylabel('$\sigma_w$','Interpreter','latex')
figPref = ['N100M20_sigW_','h',num2str(h),'_',date];
saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
print('-depsc',[saveFolder,filesep,figPref,'.eps'])

%% How sparsity and distribution of W change with different input, such as powerlaw, lognorml or exponential
dFolder = '/Users/shan/Dropbox/olfactionProject/data/powerLaw';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

fName = 'int_N100_R20_S4_sig1.5_2018-06-29.mat';
load(fullfile(dFolder,fName))

N=100;
R=20;
sp =4;
alp = 1.5;  % the exponent of power law

% select an example matrix 
inx = 4;
%==========================================================
% plot the histogram of the w
%==========================================================
temp = allMat(:,inx);

figure
set(gcf,'Renderer','painters')
histogram(temp(temp < 2),'Normalization','probability')
xlim([-80,5])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('probability')
prefix = ['gcmi_power_exampleW_hist',num2str(N),'M',num2str(R),'sp',num2str(sp),...
    '_alpha',num2str(alp),'_',date];
fileNameEps = [outFolder,filesep,prefix,'.eps'];
fileNameFig = [outFolder,filesep,prefix,'.fig'];
% print('-depsc','-painters',fileNameEps)
print(gcf,'-depsc',fileNameEps)
saveas(gcf,fileNameFig)



% zoom in of the active part
figure
set(gcf,'Renderer','painters')
% histogram(temp(temp > -10),20,'Normalization','probability')
histogram(temp(abs(temp) < 10),'Normalization','probability')
xlim([-10,1])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('probability')
prefix = ['gcmi_power_exampleW_acti',num2str(N),'M',num2str(R),'sp',num2str(sp),...
    '_alpha',num2str(alp),'_',date];
fileNameEps = [outFolder,filesep,prefix,'.eps'];
fileNameFig = [outFolder,filesep,prefix,'.fig'];
print(gcf,'-depsc',fileNameEps)
saveas(gcf,fileNameFig)



% heatmap of the examplified matrix
W = reshape(temp,R,N);
figure
set(gcf,'Renderer','painters')
imagesc(W,[-8,0])
xlabel('odorant')
ylabel('receptor')
prefix = ['gcmi_power_exampleW_heatmap',num2str(N),'M',num2str(R),'sp',num2str(sp),...
    '_alpha',num2str(alp),'_',date];
fileNameEps = [outFolder,filesep,prefix,'.pdf'];
fileNameFig = [outFolder,filesep,prefix,'.fig'];
print(gcf,'-dpdf',fileNameEps)
saveas(gcf,fileNameFig)

% =========================================================================
% Compare the distribution of power law and lognormal
% =========================================================================
% parameters for lognormal
cmin = 1;
mu = 5;
sig = 2;
pdfCut = @(x) sqrt(2/pi/sig^2)/(erfc((log(cmin)-mu)/sqrt(2)/sig))*exp(-(x-mu).^2/2/sig^2);
X = log(cmin) : 0.05:15;
Y1 = normpdf(X,mu,sig);
% Y1 = pdfCut(X);
% parameters for power law
alp = 1.2;
pd = @(X) (alp -1)*cmin^(alp -1)*exp((1-alp)*X);
Y2 = pd(X);
figure
hold on
plot(X,Y1)
plot(X,Y2)
lg = legend('lognormal,\mu = 5,\sigma = 2','power law, \alpha = 1.3','Location','northoutside');
set(lg,'FontSize',18)
legend boxoff
xlim([log(cmin),15])
ylim([0,0.2])

% yl = get(gca,'YLim');
% plot([log(cmin);log(cmin)],[yl(1);yl(2)],'k--','LineWidth',1.5)
hold off
xlabel('$\ln(c)$','Interpreter','latex')
ylabel('pdf')

z1 = trapz(X,Y2)

% =========================================================================
% Compare the distribution of power law and lognormal using gcmi simulation
% data, this is just a rough 
% =========================================================================
logNormFile = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_ref/int_N100_R20_S5_sig2_2018-06-12.mat';
powerLawFile = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_powerlaw_noReg/int_N100_R20_S5_sig1.5_2018-06-29.mat';
load(logNormFile)
Wln = allMat(allMat > -10);

load(powerLawFile)
Wpl = allMat(abs(allMat) < 10);

% histogram
figure
histogram(Wln-1,40,'DisplayStyle','stairs','LineWidth',2,'Normalization','pdf')
hold on
histogram(Wpl,40,'DisplayStyle','stairs','LineWidth',2,'Normalization','pdf')
hold off
legend('lognormal','power law')
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('pdf')

%% How differential entropy changes with the ratio of inhibition if we consider both inhibition and excitation
% the data in this section is from Shanshan's simulation
% dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_inhi';
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_inhi/N50M10S2sig2_1013'; % new data
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
allFile = dir([dFolder,filesep,'*.mat']);
allfiles = {allFile.name}';

nRecp = 10;  %need to specify this part for different data set
str_marker = '2018-10-13'; % old data, use 2018-06-22
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(allfiles,str,'once'));
files = allfiles(FIND(str_marker));

% default graphics settings
defaultGraphicsSetttings
allRatio = zeros(length(files),1);
stdRatio = zeros(length(files),1);
meanFmin = zeros(length(files),1);
stdFmin = zeros(length(files),1);
allr0 = zeros(length(files),1);
exciEle = zeros(length(files),2); %mean and std
inhiEle = zeros(length(files),2); %mean and std

% =====================================================
% make a movie of how the inhibitory and excitatory wij change
% =====================================================
% set the color
lpurple = [230,61,144]/256;
lblue = [46,167,224]/256;
myGreen = [121,191,132]/256;
myRed = [234,85,65]/256;
for i0 = 1:length(files)
    str = '(?<= *alp)[\d.]+(?=_)';
    allr0(i0) = str2num(char(regexp(files{i0},str,'match')));
    load(files{i0})
    allRatio(i0) = sum(allSign(:)  < 0)/length(allMat(:));
    stdRatio(i0) = std(sum(allSign < 0,1)/length(allMat(:,1)));
    meanFmin(i0) = - mean(allfmin(allfmin < 2));
    stdFmin(i0) = std(allfmin(allfmin < 10));
    exciEle(i0,:) =  [mean(log(allMat(allSign > 0))),std(log(allMat(allSign > 0)))];
    inhiEle(i0,:) =  [mean(log(allMat(allSign < 0))),std(log(allMat(allSign < 0)))];

    figure
    set(gcf,'Renderer','painters')
    histogram(log(allMat(allSign > 0)),'EdgeColor','none','FaceColor',lblue)
    hold on
    histogram(log(allMat(allSign  < 0)),'EdgeColor','none','FaceColor',lpurple)
    
    % add a patch
    YL = 600;
    YH = 900;
    sft = 1;  % this is due to change axis range
    xp1 = sft + [3.8,4.5,4.5,3.8];
    yp1 = [YL,YL,YH,YH];
    v1 = [sft + 3.5, YL; sft + 4.2, YL;sft + 4.2, YH;sft + 3.5, YH];
    f1 = [1 2 3 4];
    patch('Faces',f1,'Vertices',v1,...
    'EdgeColor','k','FaceColor','none','LineWidth',2);
    
    v2 = [sft + 3.52, YL;sft + 4.18, YL;sft + 4.18, YL + (YH-YL)*(allr0(i0));...
        sft + 3.52, YL + (YH-YL)*(allr0(i0))];
    f2 = [1 2 3 4];
    patch('Faces',f2,'Vertices',v2,...
    'EdgeColor','none','FaceColor',myRed);
    hold off
    xp = 2.2+sft;
    yp = YH+40;
    txt = ['$r_0$ =', num2str(allr0(i0))];
    text(xp,yp,txt,'Interpreter','latex','FontSize',20)

%     xlim([-3 5])
%     ylim([0,800])
    xlim([-3,6])
    ylim([0,1050])

    lg = legend('excitation','inhibition','Location','northwest');
    set(lg,'FontSize',16)
    legend boxoff
    xlabel('$\ln(w)$','Interpreter','latex')
    ylabel('counts')
    prefix = ['inhi_acti_histW_N50R',num2str(nRecp),'_',num2str(i0),'.png'];
    print(gcf,'-r600','-dpng',fullfile(sFolder,prefix))
%     print(gcf,'-depsc',fullfile(sFolder,prefix))
end

% ==================================================================
% make a movie to show how differential entropy changes with sparsity
%==================================================================

for i0  = 1:length(meanFmin)
    figure
    set(gcf,'Renderer','painters')
    xlim([0,1])
    ylim([-1,0])
    hold on
    errorbar(allr0(1:i0),meanFmin(1:i0),stdFmin(1:i0),'o-','MarkerSize',12,...
        'Color',Bu(10,:), 'MarkerFaceColor',Bu(10,:),'LineWidth',2,'CapSize',0)
    xlabel('basal activity')
    ylabel('differential entropy')
    
    % add a patch
    YL = -0.95;
    YH = -0.65;
    XL = 0.65;
    XH = 0.72;

    v1 = [XL,YL; XH, YL;XH, YH;XL, YH];
    f1 = [1 2 3 4];
    patch('Faces',f1,'Vertices',v1,'EdgeColor','k','FaceColor','none','LineWidth',2);
    
    v2 = [XL, YL; XH, YL; XH, YL + (YH-YL)*(allRatio(i0));XL, YL + (YH-YL)*(allRatio(i0))];
    f2 = [1 2 3 4];
    patch('Faces',f2,'Vertices',v2,'EdgeColor','none','FaceColor',myRed);
    hold off
    xp = 0.6;
    yp = YH+0.05;
    txt = ['$\rho_{\mathrm{inhi}}$ =', num2str(round(allRatio(i0)*100)/100)];
    text(xp,yp,txt,'Interpreter','latex','FontSize',20)
    prefix = ['inhi_fmin_histW_N50R',num2str(nRecp),'_',num2str(i0),'.png'];
    print(gcf,'-r600','-dpng',fullfile(sFolder,prefix))
end

% the final fmin as a function of r0
figure
errorbar(allr0',meanFmin,stdFmin,'o-','MarkerSize',12,'Color',Bu(10,:), ...
    'MarkerFaceColor',Bu(10,:),'LineWidth',2,'CapSize',0)
lg = legend('N=50,M=10,sp=2,\sigma_c = 2','Location','southeast');
set(lg,'FontSize',20)
legend boxoff
xlabel('basal activity')
ylabel('differential entropy')
prefix = ['inhi_fmin_histW_N50R',num2str(nRecp),];
print(gcf,'-depsc',fullfile(sFolder,[prefix,'.eps']))
saveas(gcf,fullfile(sFolder,[prefix,'.fig']))

% ==================================================================
% plot the inhibitory ratio vs r0
% ==================================================================
figure
hold on
errorbar(allr0',allRatio,stdRatio,'o-','MarkerSize',12,'Color',Bu(10,:), ...
    'MarkerFaceColor',Bu(10,:),'LineWidth',2,'CapSize',0)
lg = legend('N=50,M=10,sp=2,\sigma_c = 2','Location','northwest');
set(lg,'FontSize',20)
legend boxoff
xlim([0,1]);ylim([0,1]);
plot([0;1],[0;1],'k--','LineWidth',1.5)
hold off
box on
xlabel('basal activity')
ylabel('inhibitory ratio')
prefix = ['inhi_inhiRatio_histW_N50R',num2str(nRecp),];
print(gcf,'-depsc',fullfile(sFolder,[prefix,'.eps']))
saveas(gcf,fullfile(sFolder,[prefix,'.fig']))
% ==================================================================
% heatmap of an example w
%==================================================================
selFile = 'int_N50_R10_S2_sig2_alp0.5_2018-06-22.mat';
load(fullfile(dFolder,selFile));
w = reshape(allMat(:,1),[10,50]);
w(allSign(:,1)<0) = nan;

figure
imagesc(log(w),[-2,2])
set(gca,'FontSize',16)
colormap
prefix = 'inhi_exampW_N50M10_sp2_heatmap_active';
print('-depsc',fullfile(sFolder,[prefix,'.eps']));
saveas(gcf,fullfile(sFolder,[prefix,'.fig']))

% ==================================================
%  plot three examples of W distribution 

% =====================================================


sInx = [2,6,10];
sr0  = [0.1,0.5,0.8];
for i0 = 1:length(sInx)
    load(files{sInx(i0)})
    figure
    hold on
    histogram(log(allMat(allSign > 0)),40,'EdgeColor','none','FaceColor',lpurple)
    histogram(log(allMat(allSign  < 0)),40,'EdgeColor','none','FaceColor',lblue)
    hold off
    xlim([-4,3])
    legend('excitation','inhibition')
    legend boxoff
    xlabel('$\ln(w)$','Interpreter','latex')
    ylabel('counts')
    
%     fileNameEps = ['N50_R20_sp5_sig2_inhi_wdist_r0_',num2str(sr0(i0)),...
%         '_',date,'.eps'];
%     fileNameFig = ['N50_R20_sp5_sig2_inhi_wdist_r0_',num2str(sr0(i0)),...
%         '_',date,'.fig'];
%     print('-depsc',fileNameEps)
%     saveas(gcf,fileNameFig)
end

% ============ plot an example matrix, use r0 =0.5 ===============
% ================================================================
sInx = 6;
load(files{sInx})
w = reshape(allMat(:,1),[10,50]);
%only consider the active part
w1 = log(w);
w1(allSign(:,1) < 0) = nan;  %get rid of inhibiotry response
hh = heatmap(w1);
hh.ColorLimits = [-5,2];
w1(w1) = log(w(allSign(:,1) > 0));
figure


imagesc(w,[-3,2])
figure
for i0 = 1:length(files)
    str = '(?<= *alp)[\d.]+(?=_)';
    allr0(i0) = str2num(char(regexp(files{i0},str,'match')));
    load(files{i0})
    allRatio(i0) = sum(allSign(:)  < 0)/length(allMat(:));
    stdRatio(i0) = std(sum(allSign < 0,1)/length(allMat(:,1)));
    meanFmin(i0) = - mean(allfmin(allfmin < 2));
    stdFmin(i0) = std(allfmin(allfmin < 2));
    exciEle(i0,:) =  [mean(log(allMat(allSign > 0))),std(log(allMat(allSign > 0)))];
    inhiEle(i0,:) =  [mean(log(allMat(allSign < 0))),std(log(allMat(allSign < 0)))];
%     histogram(log(allMat(allSign > 0)),'Normalization','probability')
%     histogram(log(allMat(allSign  < 0)),'Normalization','probability')
    subplot(4,4,i0)
    histogram(log(allMat(allSign > 0)),40)
    hold on
    histogram(log(allMat(allSign  < 0)),40)
    hold off
%     xlim([-5 3])
    if i0==1
        legend('excitation','inhibition')
        legend boxoff
    end
    xlabel('$\ln(w)$','Interpreter','latex')
    ylabel('counts')
end

% ============================================================
% plot how fmin changes with spontaneous activity
% ============================================================
[uniqr0,inx] = sort(allr0);
figure
errorbar(uniqr0',meanFmin(inx),stdFmin(inx),'o-','MarkerSize',12,'LineWidth',...
    2,'Color',Bu(10,:),'MarkerFaceColor',Bu(10,:))
lg = legend('N=50,M=10,sp=2,\sigma_c = 2');
set(lg,'FontSize',20)
legend boxoff
set(gca,'XTick',0:0.2:1)
xlabel('$r_0$','Interpreter','latex')
% ylabel('$f_{min}$','Interpreter','latex')
ylabel('differential entropy')

fileNameEps = ['N50_R',num2str(nRecp),'_sp5_sig2_inhi_fmin_',date,'.eps'];
fileNameFig = ['N50_R',num2str(nRecp),'_sp5_sig2_inhi_fmin_',date,'.fig'];
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)

% ============================================================
% how inhibitory ratio changes with r_0
% ============================================================
figure
hold on
errorbar(uniqr0',allRatio(inx),stdRatio(inx),'o-','MarkerSize',12,'LineWidth',...
    2,'Color',Bu(10,:),'MarkerFaceColor',Bu(10,:))
plot([0,1],[0,1],'k--','LineWidth',1)
hold off
lg = legend('N=50,M=10,sp=2,\sigma_c = 2');
set(lg,'FontSize',20)
legend boxoff
box on
set(gca,'XTick',0:0.2:1)
xlabel('$r_0$','Interpreter','latex')
% ylabel('$\rho$','Interpreter','latex')
ylabel('inhibitory ratio')

fileNameEps = ['N50_R',num2str(nRecp),'sp5_sig2_inhi_ratio_',date,'.eps'];
fileNameFig = ['N50_R',num2str(nRecp),'sp5_sig2_inhi_ratio_',date,'.fig'];
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)


% distribution of excitatory and inhibitory elements
% firs the mean value
figure
hold on
plot(uniqr0',exciEle(inx,1),'o-','MarkerSize',12,'LineWidth',2)
plot(uniqr0',inhiEle(inx,1),'o-','MarkerSize',12,'LineWidth',2)
hold off
legend('exci','inhi')
legend boxoff
xlabel('$r_0$','Interpreter','latex')
ylabel('$\mu_w$','Interpreter','latex')
fileNameEps = ['N50_R20_sp5_sig2_inhi_meanVal_',date,'.eps'];
fileNameFig = ['N50_R20_sp5_sig2_inhi_meanVal_',date,'.fig'];
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)

figure
hold on
plot(uniqr0',exciEle(inx,2),'o-','MarkerSize',12,'LineWidth',2)
plot(uniqr0',inhiEle(inx,2),'o-','MarkerSize',12,'LineWidth',2)
hold off
legend('exci','inhi')
legend boxoff
ylim([0,2])
xlabel('$r_0$','Interpreter','latex')
ylabel('$\sigma_w$','Interpreter','latex')
fileNameEps = ['N50_R20_sp5_sig2_inhi_sigW_',date,'.eps'];
fileNameFig = ['N50_R20_sp5_sig2_inhi_sigW_',date,'.fig'];
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)


% di
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

%% This part is for the power law distribution of input
% dataFolder = '/Users/shan/Dropbox/olfactionProject/data/';
dataFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_powerlaw_diffBeta';
figFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
fileNameCell = {allFile.name}';

num = 40;   % repeat of simulation
% this set of data happened to with the same alpha
% str_marker = '(?<= *_sig)[2]+(?=_)';   %folder with this string contains the data we need
str_marker = 'N50_R13';   %folder with this string contains the data we need
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fileNameCell,str,'once'));
files = fileNameCell(FIND(str_marker));

%============
% threshold to register the active elements
%======
thd = -10;


s1 = '(?<= *N)[\d.]+(?=_)';
s2 = '(?<= *_R)[\d]+(?=_)';
s3 = '(?<= *_S)[\d]+(?=_)';
s5 = '(?<= *sig)[\d.]+(?=_)';

numOdor = zeros(length(files),1);
numRecp = zeros(length(files),1);
spInput = zeros(length(files),1);
sig = zeros(length(files),1);
spW = zeros(length(files),1);
% spWMat = zeros(length(files), num); % store the sparsity of individual matrix
stdSpW = zeros(length(files),1);
meanFmin = zeros(length(files),1);
stdFmin = zeros(length(files),1);
% figure
for i0 = 1:length(files)
    load(fullfile(dataFolder,files{i0}))
    numOdor(i0) = str2num(char(regexp(files{i0},s1,'match')));
    numRecp(i0) = str2num(char(regexp(files{i0},s2,'match')));
    spInput(i0) = str2num(char(regexp(files{i0},s3,'match')));
    sig(i0) = str2num(char(regexp(files{i0},s5,'match')));
    
    t = allMat(:);
    spW(i0) = sum(t>thd)/length(t);
    stdSpW(i0) = std(sum(allMat > thd,1)/size(allMat,1));
    
    meanFmin(i0) = mean(-allfmin); 
    stdFmin(i0) = std(-allfmin);
end

% plot
figure
[~,ix] = sort(round(sig*10));
errorbar(sig(ix),spW(ix),stdSpW(ix),'o-','MarkerSize',12,'MarkerFaceColor',...
    myBlue,'Color',myBlue,'CapSize',0,'LineWidth',2)
xlabel('$\beta$','Interpreter','latex')
ylabel('$\rho_w$','Interpreter','latex')
ylim([0.55,0.83])
prefix = ['powerLaw_spW_diffBeta_',date];
fileNameEps = [figFolder,filesep,prefix,'.eps'];
fileNameFig = [figFolder,filesep,prefix,'.fig'];
print('-depsc','-painters',fileNameEps)
saveas(gcf,fileNameFig)


% group the data
uniqOdor = sort(uniq(numOdor));
uniqSig = [1.2,1.5,2,2.5,3];
% uniqSig = [1.5,3];
uniqSp = sort(uniq(spInput));
sparsity = cell(length(uniqOdor),1);
spOut = cell(length(uniqSig),1);
stdSpOut = cell(length(uniqSig),1);
allAveFmin = cell(length(uniqSig),1);
allStdFmin = cell(length(uniqSig),1);

for i0 = 1:length(uniqSig)
    inx = sig == uniqSig(i0);
    
    sparsity{i0} = spInput(inx);
    spOut{i0} = spW(inx);
    stdSpOut{i0} = stdSpW(inx);
    allAveFmin{i0} = meanFmin(inx);
    allStdFmin{i0} = stdFmin(inx);
end

%plot the output sparsity and input sparisty
defaultGraphicsSetttings
figure
hold on
for i0 = 1:length(uniqSig)
    plot(sparsity{i0},spOut{i0},'o-','MarkerSize',12,'Color',Bu(4 + i0,:),...
    'MarkerFaceColor',Bu(4 + i0,:),'LineWidth',2)
end
hold off
legend(['\beta = ',num2str(uniqSig(1))],['\beta = ',num2str(uniqSig(2))])
legend boxoff

% errorbar(sparsity{2},spOut{2},stdSpOut{2},'o-','MarkerSize',12,'Color',Bu(10,:),...
%     'MarkerFaceColor',Bu(10,:),'LineWidth',2,'CapSize',0)
% legend('\alpha = 2')
ylim([0.3,0.85])
box on
legend boxoff
xlabel('input sparsity')
ylabel('sparsity of w')

% fileNameEps = [figFolder,filesep,'powerLaw_spW_spc_diff_alp_gcmi.eps'];
% fileNameFig = [figFolder,filesep,'powerLaw_spW_spc_diff_alp_gcmi.fig'];
% prefix = ['powerLaw_spW_spc_alp',num2str(uniqSig(2)),'_gcmi'];
prefix = ['powerLaw_spW_spc_diffAlp_gcmi',date];

fileNameEps = [figFolder,filesep,prefix,'.eps'];
fileNameFig = [figFolder,filesep,prefix,'.fig'];
print('-depsc','-painters',fileNameEps)
saveas(gcf,fileNameFig)


% ================================================================
% plot how differential entropy change with sp for different alpha
% ================================================================
figure
hold on
for i0 = 1:length(uniqSig)
    errorbar(sparsity{i0},allAveFmin{i0},allStdFmin{i0},'o-','MarkerSize',12,...
        'Color',Bu(4 + 3*i0,:),'MarkerFaceColor',Bu(4 + 3*i0,:),'LineWidth',2)
end
hold off
legend(['\beta = ',num2str(uniqSig(1))],['\beta = ',num2str(uniqSig(2))])
legend boxoff
box on
% ylim([-30,-15])

xlabel('input sparsity')
ylabel('differential entropy')

prefix = ['expon_fmin_spc_diffAlp_gcmi',date];
fileNameEps = [figFolder,filesep,prefix,'.eps'];
fileNameFig = [figFolder,filesep,prefix,'.fig'];
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)

% ===============================================================
% heatmap and histogram of an example 
% ===============================================================
selFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_powerlaw';
selFile = 'Gcmi_power_N100_R20_S5_sig1.5_2018-10-15.mat';
load(fullfile(selFolder,selFile))

N = 100;
M = 20;
sp = 5;
alp = 1.5;

ix = 9;  % select on of the matrix
w = reshape(allMat(:,ix),[M,N]);

% heatmap
figure
set(gcf,'renderer','Painters')
imagesc(w,[-6,3])
set(gca,'TickLength',[0,0])
colormap(jet);
c = colorbar; 
xlabel('Odorant')
ylabel('Receptor')
prefix = ['powerlaw_exampleW_N',num2str(N),'M',num2str(M),'alp',num2str(alp),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[figFolder,filesep,prefix,'.fig'])
print('-depsc',[figFolder,filesep,prefix,'.eps'])


% histogram
figure
hold on
set(gcf,'renderer','Painters')
h1 = histogram(w(:),60,'Normalization','probability');
h1.FaceColor = lBu; h1.FaceAlpha = 0.4; h1.EdgeColor = 'none';
stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
    [0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',2)

box on
hold off
set(gca,'Layer','top')
xlim([-80,5])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('probability')
prefix = ['power_exampleW_N',num2str(N),'M',num2str(M),'alp',num2str(alp),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[figFolder,filesep,prefix,'.fig'])
print('-depsc',[figFolder,filesep,prefix,'.eps'])


% enlarged
figure
set(gcf,'renderer','Painters')
hold on
h2 = histogram(w(abs(w)<8),20,'Normalization','pdf');
h2.FaceColor = lBu; h2.FaceAlpha = 0.4; h2.EdgeColor = 'none';
stairs([h2.BinEdges(1),h2.BinEdges,h2.BinEdges(end)],...
    [0,h2.Values,h2.Values(end),0],'Color',dpBu,'LineWidth',2)

% get the normal fit parameter
% pd = fitdist(w(abs(w)<6),'normal');
% X = -6:0.05:6;
% Y = normpdf(X,pd.mean,pd.sigma);
% plot(X,Y,'Color',Or,'LineWidth',4)
hold off
box on
% legend('active w','Gassian fit')
% legend boxoff
set(gca,'XLim',[-6,2])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
prefix = ['power_exampleActiveW_N',num2str(N),'M',num2str(M),'alp',num2str(alp),....
    'sp',num2str(sp),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])



% =========================================================
% plot the tail distribution of w elements
% =========================================================
% select an examplified data
selFile = 'Gcmi_power_N100_R20_S5_sig1.5_2018-10-10.mat';
N =100;
M = 20;
sp = 5;
alp = 1.5;
load(fullfile(dataFolder,selFile))

% plot the histogram of Wij
figure
temp = allMat(:);
histogram(temp(temp < 5),'Normalization','probability')
xlim([-80,5])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('probability')
prefix = ['gcmi_power_exampleW_N',num2str(N),'M',num2str(M),'sp',num2str(sp),...
    '_alpha',num2str(alp),'_',date];
fileNameEps = [figFolder,filesep,prefix,'.eps'];
fileNameFig = [figFolder,filesep,prefix,'.fig'];
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)


% zoom in on the active part
figure
temp = allMat(:);
histogram(temp(abs(temp) < 10),'Normalization','probability')
xlim([-10,0])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('probability')
prefix = ['gcmi_power_exampleW_active_N',num2str(N),'M',num2str(M),'sp',num2str(sp),...
    '_alpha',num2str(alp),'_',date];
fileNameEps = [figFolder,filesep,prefix,'.eps'];
fileNameFig = [figFolder,filesep,prefix,'.fig'];
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)

% test if the left tail is power law distributed
Wact = temp(abs(temp) < 9);
[F,X] = ecdf(Wact);
inx = 1:20:length(X);
plot(X(inx),log10(F(inx)))
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('cdf')
% plot the heatmap of an example matrix
[~,inx] = sort(allfmin,'ascend');
w = reshape(allMat(:,inx(1)),[M,N]);
figure
imagesc(w,[-8,1])
set(gca,'FontSize',16)
xlabel('odorant')
ylabel('receptor')
prefix = ['gcmi_power_exampleW_heatMap',num2str(N),'M',num2str(M),'sp',num2str(sp),...
    '_alpha',num2str(alp),'_',date];
fileNameEps = [figFolder,filesep,prefix,'.eps'];
fileNameFig = [figFolder,filesep,prefix,'.fig'];
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)






% ===================   plot how results depend on alpha ============

str_marker = '2018-06-13';       %folder with this string contains the data we need
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

%% This part show the deviation from Gauss of optimal W when odor is power law distribution
% here we compare the skewness and kurtosis and compare with corresponding
% Gausssian distribtuion with finite size (same mean and standard deviaiton)

selFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_powerlaw10102018';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

selFile1 = 'Gcmi_power_N100_R20_S3_sig1.5_2018-10-10.mat';
load(fullfile(selFolder,selFile1))
bt = 1.5;   %exponent of power law distribution

L = length(allfmin);
allKt = zeros(L,1);
allSk = zeros(L,1);
GaussCmp = zeros(L,2); %store the Gaussian comparision
thd = -8; % the choice of this threshold should be carefull
for i0 = 1:L
    temp = allMat(:,i0);
    sensiW = temp(temp > thd & temp < 5);
    allSk(i0) = skewness(sensiW);
    allKt(i0) = kurtosis(sensiW);  
    
    % generate Gaussian data
    gd = std(sensiW)*randn(length(sensiW),1) + mean(sensiW);
    GaussCmp(i0,1) = skewness(gd);
    GaussCmp(i0,2) = kurtosis(gd);
end

% export data for plot in R
dataFileName = fullfile(sFolder,['SkewKurt_powerlaw_N100M20S5_beta',num2str(bt),'.mat']);
save(dataFileName,'GaussCmp','allSk','allKt','bt')

% scatter plot to compare with Gaussian distribution
figure
hold on
scatter(GaussCmp(:,1),allSk,40,'filled')
xlim([-0.3,0.3])
ylim([-0.3,0.3])
plot([-0.3;0.3],[-0.3,0.3],'k--')
box on
xlabel('Gaussian')
ylabel('simulation')

% kurtosis
figure
hold on
scatter(GaussCmp(:,2)-3,allKt-3,40,'filled')
xlim([-1.1,0.4])
ylim([-1.1,0.4])
plot([-1.1,0.4],[-1.1,0.4],'k--')
box on
xlabel('Gaussian')
ylabel('simulation')


% barplot with erorr bar, comparison
allError= zeros(2,2);
summMean = zeros(2,2);

allError(1,:) = std(GaussCmp,0,1);
allError(2,:) = [std(allSk),std(allKt)];
summMean(1,:) = [mean(GaussCmp(:,1)),mean(GaussCmp(:,2)-3)];
summMean(2,:) = [mean(allSk),mean(allKt)-3];

figure
h = barwitherr(allError', summMean');% Plot with errorbars

set(gca,'XTickLabel',{'Skewness','Excess Kurtosis'})
legend('Gaussian','Simulation')
ylabel('Value')
set(h(1),'FaceColor','k');


% ================================================================
% show how the skewness change with beta of power law distribution
% ================================================================
dataFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_powerlaw_diffBeta20190128';
figFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
fileNameCell = {allFile.name}';

num = 40;   % repeat of simulation
% this set of data happened to with the same alpha
% str_marker = '(?<= *_sig)[2]+(?=_)';   %folder with this string contains the data we need
str_marker = 'N50_R13';   %folder with this string contains the data we need
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fileNameCell,str,'once'));
files = fileNameCell(FIND(str_marker));

% threshold to register the active elements
thd = -8;

s1 = '(?<= *N)[\d.]+(?=_)';
s2 = '(?<= *_R)[\d]+(?=_)';
s3 = '(?<= *_S)[\d]+(?=_)';
s5 = '(?<= *sig)[\d.]+(?=_)';


GaussSk = zeros(length(files),num);
GaussKt = zeros(length(files),num);
allSk = zeros(length(files),num);
allKt = zeros(length(files),num);
fmin = zeros(length(files),2);
allSparsity = zeros(length(files),num);
for i0 = 1:length(files)
    load(fullfile(dataFolder,files{i0}))
%     numOdor(i0) = str2num(char(regexp(files{i0},s1,'match')));
%     numRecp(i0) = str2num(char(regexp(files{i0},s2,'match')));
%     spInput(i0) = str2num(char(regexp(files{i0},s3,'match')));
    sig(i0) = str2num(char(regexp(files{i0},s5,'match')));
    
%     t = allMat(:);
%     spW(i0) = sum(t>thd)/length(t);
%     stdSpW(i0) = std(sum(allMat > thd,1)/size(allMat,1));
%     
%     meanFmin(i0) = mean(-allfmin); 
%     stdFmin(i0) = std(-allfmin);
    fmin(i0,:) = [mean(-allfmin),std(-allfmin)];
    for j0 = 1:size(allMat,2)
        temp = allMat(:,j0);
        sensiW = temp(temp > thd & temp < 5);
        allSk(i0,j0) = skewness(sensiW);
        allKt(i0,j0) = kurtosis(sensiW)-3;  
        allSparsity(i0,j0) = sum(temp > thd)/length(temp(:));
    % generate Gaussian data
        gd = std(sensiW)*randn(length(sensiW),1) + mean(sensiW);
        GaussSk(i0,j0) = skewness(gd);
        GaussKt(i0,j0) = kurtosis(gd)-3;
    end
end

%plot the summary of the results
simSk = [mean(allSk,2),std(allSk,0,2)];
gaussSk = [mean(GaussSk,2),std(GaussSk,0,2)];
simKt = [mean(allKt,2),std(allKt,0,2)];
GaussKt = [mean(GaussKt,2),std(GaussKt,0,2)];
spW = [mean(allSparsity,2),std(allSparsity,0,2)];
% plot how skewness change with beta
figure
hold on
[~,inx] = sort(round(sig*100));
errorbar(sig(inx),simSk(inx,1),simSk(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
errorbar(sig(inx),gaussSk(inx,1),gaussSk(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Gr(7,:),'Color',Gr(7,:),'LineWidth',2,'CapSize',0)
lg = legend('simulation','Gaussian');
set(lg,'FontSize',20)
ylim([-0.5,0.5])
xlabel('$\beta$','interpreter','latex')
ylabel('Skewness')
box on
prefix = ['Skewness_powerlaw_beta_N50M13_',date];
saveas(gcf,[figFolder,filesep,prefix,'.fig'])
print('-depsc',[figFolder,filesep,prefix,'.eps'])


% plot how kurtosis change with beta
figure
hold on
[~,inx] = sort(round(sig*100));
errorbar(sig(inx),simKt(inx,1),simKt(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
errorbar(sig(inx),GaussKt(inx,1),GaussKt(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Gr(7,:),'Color',Gr(7,:),'LineWidth',2,'CapSize',0)
lg = legend('simulation','Gaussian');
set(lg,'FontSize',20)
box on
ylim([-1.5,1])
xlabel('$\beta$','interpreter','latex')
ylabel('Excess kurtosis')
prefix = ['Kurtosis_powerlaw_beta_N50M13_',date];
saveas(gcf,[figFolder,filesep,prefix,'.fig'])
print('-depsc',[figFolder,filesep,prefix,'.eps'])


% plot how sparsity change with alpha
figure
hold on
[~,inx] = sort(round(sig*100));
errorbar(sig(inx),spW(inx,1),spW(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
box on
ylim([0.1,0.9])
xlabel('$\beta$','interpreter','latex')
ylabel('$\rho_w$','interpreter','latex')
prefix = ['diffEntr_skew_alpha_N50M13_',date];
saveas(gcf,[figFolder,filesep,prefix,'.fig'])
print('-depsc',[figFolder,filesep,prefix,'.eps'])


% plot how optimal differential entropy change with alpha
figure
hold on
[~,inx] = sort(round(sig*100));
errorbar(sig(inx),fmin(inx,1),fmin(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
box on
ylim([-15,0])
xlabel('$\beta$','interpreter','latex')
ylabel('differential entropy')
prefix = ['diffEntr_skew_alpha_N50M13_',date];
saveas(gcf,[figFolder,filesep,prefix,'.fig'])
print('-depsc',[figFolder,filesep,prefix,'.eps'])

% ========================================================
% plot example histogram
% ===================n=====================================
bt = 1.2;
figure
histogram(sensiW,'Normalization','pdf')
legend(['\beta =',num2str(bt)])
legend boxoff
xlabel('$\ln(w)$','interpreter','latex')
ylabel('pdf')
xlim([-8,-2])

%% Show the skewness of W change with input skewness
% here we use skewed Gaussian distribution as the input

dataFolder = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/code/data/gcmi_skew/N100R20sig3_noReg';
figFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
fileNameCell = {allFile.name}';

num = 40;   % repeat of simulation
% this set of data happened to with the same alpha
% str_marker = '(?<= *_sig)[2]+(?=_)';   %folder with this string contains the data we need
str_marker = '2019-02-04';   %folder with this string contains the data we need
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fileNameCell,str,'once'));
files = fileNameCell(FIND(str_marker));

N = 100; 
M = 20; 
sig = 3;  


% threshold to register the active elements
thd = -5;
UB = 5;  % get rid off ultrasensitivity

s1 = '(?<= *N)[\d.]+(?=_)';
s2 = '(?<= *_R)[\d]+(?=_)';
s3 = '(?<= *_S)[\d]+(?=_)';
s5 = '(?<= *sig)[\d.]+(?=_)';
s6 = '(?<= *alp)[-,\d]+(?=_)';


GaussSk = nan(length(files),num);
GaussKt = nan(length(files),num);
allSk = nan(length(files),num);
allKt = nan(length(files),num);
fmin = nan(length(files),2);
allSparsity = nan(length(files),num);
allCDF = cell(length(files),1);  % store the emprical cdf
% X = -5:0.005:10;


alp = zeros(length(files),1);  % store the skewness parameters
for i0 = 1:length(files)
    load(fullfile(dataFolder,files{i0}))

    alp(i0) = str2num(char(regexp(files{i0},s6,'match')));
    
%     t = allMat(:);
%     spW(i0) = sum(t>thd)/length(t);
%     stdSpW(i0) = std(sum(allMat > thd,1)/size(allMat,1));
%     
%     meanFmin(i0) = mean(-allfmin); 
%     stdFmin(i0) = std(-allfmin);
%     storeCDF = zeros(length(X),size(allMat,2));
    fmin(i0,:) = [mean(-allfmin),std(-allfmin)];
    [f,X] = ecdf(allMat(allMat > thd));
    allCDF{i0} = [f,X];
    for j0 = 1:size(allMat,2)
        temp = allMat(:,j0);
        allSparsity(i0,j0) = sum(temp > thd)/length(temp(:));
        
        sensiW = temp(temp > thd & temp <=UB);  % here the upper limit should be checked
%         [f, X] = ecdf(temp(temp > thd));  % we can also use only part of the all the sensitive 
        allSk(i0,j0) = skewness(sensiW);
        allKt(i0,j0) = kurtosis(sensiW)-3;  
    
    % generate Gaussian data
        gd = std(sensiW)*randn(length(sensiW),1) + mean(sensiW);
        GaussSk(i0,j0) = skewness(gd);
        GaussKt(i0,j0) = kurtosis(gd)-3;
    end
end

%plot the summary of the results
simSk = [nanmean(allSk,2),nanstd(allSk,0,2)];
gaussSk = [nanmean(GaussSk,2),nanstd(GaussSk,0,2)];
simKt = [nanmean(allKt,2),nanstd(allKt,0,2)];
GaussKt = [nanmean(GaussKt,2),nanstd(GaussKt,0,2)];
spW = [nanmean(allSparsity,2),nanstd(allSparsity,0,2)];

% plot how skewness change with beta
figure
hold on
[~,inx] = sort(round(alp*100));
errorbar(alp(inx),simSk(inx,1),simSk(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
errorbar(alp(inx),gaussSk(inx,1),gaussSk(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Gr(7,:),'Color',Gr(7,:),'LineWidth',2,'CapSize',0)
lg = legend('simulation','Gaussian');
set(lg,'FontSize',20)
ylim([-0.7,1])
xlabel('$\alpha$','interpreter','latex')
ylabel('Skewness')
box on
prefix = ['Skewness_SkewGauss_alp_N50M13_',date];
saveas(gcf,[figFolder,filesep,prefix,'.fig'])
print('-depsc',[figFolder,filesep,prefix,'.eps'])


% plot how kurtosis change with beta
figure
hold on
[~,inx] = sort(round(alp*100));
errorbar(alp(inx),simKt(inx,1),simKt(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
errorbar(alp(inx),GaussKt(inx,1),GaussKt(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Gr(7,:),'Color',Gr(7,:),'LineWidth',2,'CapSize',0)
lg = legend('simulation','Gaussian');
set(lg,'FontSize',20)
box on
ylim([-1,1.5])
xlabel('$\alpha$','interpreter','latex')
ylabel('Excess kurtosis')
prefix = ['Kurtosis_powerlaw_beta_N50M13_',date];
saveas(gcf,[figFolder,filesep,prefix,'.fig'])
print('-depsc',[figFolder,filesep,prefix,'.eps'])

% plot how sparsity change with alpha
figure
hold on
[~,inx] = sort(round(alp*100));
errorbar(alp(inx),spW(inx,1),spW(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
box on
ylim([0.4,0.7])
xlabel('$\alpha$','interpreter','latex')
ylabel('$\rho_w$','interpreter','latex')
prefix = ['sparsity_skew_alpha_N',num2str(N),'M',num2str(M),'sig',num2str(sig),'_',date];

saveas(gcf,[figFolder,filesep,prefix,'.fig'])
print('-depsc',[figFolder,filesep,prefix,'.eps'])


% plot how optimal differential entropy change with alpha
figure
hold on
[~,inx] = sort(round(alp*100));
errorbar(alp(inx),fmin(inx,1),fmin(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
box on
ylim([-9,-5])
xlabel('$\alpha$','interpreter','latex')
ylabel('differential entropy')
prefix = ['diffEntr_skew_alpha_N',num2str(N),'M',num2str(M),'sig',num2str(sig),'_',date];
saveas(gcf,[figFolder,filesep,prefix,'.fig'])
print('-depsc',[figFolder,filesep,prefix,'.eps'])

% ========================================
% Plot the histogram of exmaple distribution
% =========================================
figure
temp = allMat(:);
thd = -5;
alp = -4;
histogram(temp(temp> thd & temp < 10),'Normalization','pdf')
legend(['\alpha = ',num2str(alp)])
legend boxoff
xlabel('$\ln(w)$','interpreter','latex')
ylabel('pdf')
xlim([thd,2])


% ===============================
% compare the empirical cdf
% ===============================
[orderAlp,inx] = sort(round(alp*100));
selAlp =  [-4,0,4];  % only plot these three
figure
hold on
% for i0 = 1:length(alp)
%     plot(allCDF{inx(i0)}(:,2),allCDF{inx(i0)}(:,1))
% end
for i0 = 1:length(selAlp)
    ix = find(alp == selAlp(i0));
    % shift the position according to medium value
    mIdx = find(allCDF{ix}(:,1) >= 0.5,1,'first');
    mediamW = allCDF{ix}(mIdx,2);
    plot(allCDF{ix}(:,2) - mediamW,allCDF{ix}(:,1))
end
box on
xlim([-4,3])
lg = legend('\alpha = -4', '\alpha = 0', '\alpha = 4','Location','northwest');
set(lg,'FontSize',16)
legend boxoff

hold off
ylabel('emperical cdf')
xlabel('$\ln(w)$','interpreter','latex')


% ===============================
% K-L 
% ===============================

%% This part is for the exponential distribution of input
% dataFolder = '/Users/shan/Dropbox/olfactionProject/data/powerLaw';
dataFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_exp';
figFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
fileNameCell = {allFile.name}';


% this set of data happened to with the same alpha
% str_marker = '(?<= *_sig)[2]+(?=_)';   %folder with this string contains the data we need
str_marker = 'N100_R20';   %folder with this string contains the data we need
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fileNameCell,str,'once'));
files = fileNameCell(FIND(str_marker));

%============
% threshold to register the active elements
%======
thd = -10;


s1 = '(?<= *N)[\d.]+(?=_)';
s2 = '(?<= *_R)[\d]+(?=_)';
s3 = '(?<= *_S)[\d]+(?=_)';
s5 = '(?<= *sig)[\d.]+(?=_)';

numOdor = zeros(length(files),1);
numRecp = zeros(length(files),1);
spInput = zeros(length(files),1);
sig = zeros(length(files),1);
spW = zeros(length(files),1);
stdSpW = zeros(length(files),1);
meanFmin = zeros(length(files),1);
stdFmin = zeros(length(files),1);
% figure
for i0 = 1:length(files)
    load(fullfile(dataFolder,files{i0}))
    numOdor(i0) = str2num(char(regexp(files{i0},s1,'match')));
    numRecp(i0) = str2num(char(regexp(files{i0},s2,'match')));
    spInput(i0) = str2num(char(regexp(files{i0},s3,'match')));
    sig(i0) = str2num(char(regexp(files{i0},s5,'match')));
    
    t = allMat(:);
    spW(i0) = sum(t>thd)/length(t);
    stdSpW(i0) = std(sum(allMat > thd,1)/size(allMat,1));
    
    meanFmin(i0) = mean(-allfmin); 
    stdFmin(i0) = std(-allfmin);
end

% group the data
uniqOdor = sort(uniq(numOdor));
% uniqSig = [1.5,2,2.5];
uniqSig = 5;
uniqSp = sort(uniq(spInput));
sparsity = cell(length(uniqOdor),1);
spOut = cell(length(uniqSig),1);
stdSpOut = cell(length(uniqSig),1);
allAveFmin = cell(length(uniqSig),1);
allStdFmin = cell(length(uniqSig),1);

for i0 = 1:length(uniqSig)
    inx = sig == uniqSig(i0);
    
    sparsity{i0} = spInput(inx);
    spOut{i0} = spW(inx);
    stdSpOut{i0} = stdSpW(inx);
    allAveFmin{i0} = meanFmin(inx);
    allStdFmin{i0} = stdFmin(inx);
end

%plot the output sparsity and input sparisty
defaultGraphicsSetttings
figure
hold on
for i0 = 1:length(uniqSig)
    plot(sparsity{i0},spOut{i0},'o-','MarkerSize',12,'Color',Bu(10,:),...
    'MarkerFaceColor',Bu(10,:),'LineWidth',2)
end
hold off
% legend('\alpha = 1.5','\alpha = 3')
legend('\lambda = 0.2')

% errorbar(sparsity{2},spOut{2},stdSpOut{2},'o-','MarkerSize',12,'Color',Bu(10,:),...
%     'MarkerFaceColor',Bu(10,:),'LineWidth',2,'CapSize',0)
% legend('\alpha = 2')
legend boxoff
ylim([0.3 0.55])
xlabel('input sparsity')
ylabel('sparsity of w')

% fileNameEps = [figFolder,filesep,'powerLaw_spW_spc_diff_alp_gcmi.eps'];
% fileNameFig = [figFolder,filesep,'powerLaw_spW_spc_diff_alp_gcmi.fig'];
% prefix = ['powerLaw_spW_spc_alp',num2str(uniqSig(2)),'_gcmi'];
prefix = ['expon_spW_spc_diffAlp_gcmi',date];

fileNameEps = [figFolder,filesep,prefix,'.eps'];
fileNameFig = [figFolder,filesep,prefix,'.fig'];
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)


% ================================================================
% plot how differential entropy change with sp for different alpha
% ================================================================



% ===================================================
% example data
% ===================================================
% selFile = 'int_N20_R9_S3_sig5_2018-06-29.mat';
selFile = 'int_N100_R20_S3_sig5_2018-06-29.mat';
load(fullfile(dataFolder,selFile))

% over all histogram
figure
histogram(allMat(:),'Normalization','probability')
xlim([-70,5])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('probability')

% zoom in the active part
figure
t = allMat(:);
histogram(t(t>-10),'Normalization','probability')
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('probability')

%heatmap

%% Analyze the optimal W for situation with both excitation and inhibition
% in this section, data is from Qianyi's simulation
% dFolder = '/Users/shan/Dropbox/olfactionProject/data/GcmiInhibit-06-13';
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_inhi/N50M10S2sig2_1013';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir([dFolder,filesep,'*.mat']);
files = {allFile.name}';

N = 50;
M = 10;
sig = 10;
sp = 2;

R0 = [0.05,0.08,0.1,0.12,0.14,0.16,0.18,0.25,0.3,0.4,0.5,0.6,0.7,0.82,0.84,0.88,0.9,0.95];  % the basal activity

nOdor = zeros(length(files),1);
meanRatio = nan(length(R0),length(files));
stdRatio = nan(length(R0),length(files));
meanFmin = nan(length(R0),length(files));
stdFmin = nan(length(R0),length(files));
for i0 = 1:length(files)
    str = '(?<= *basal)[\d.]+(?=_)';
%     nOdor(i0) = str2num(char(regexp(files{i0},str,'match')));
    load(fullfile(dFolder,files{i0}));
    if i0 == 1
       jstart = 0;
    else
        jstart = 1;
    end
        
    for j0 = 1:length(allMat)
        meanRatio(j0+jstart,i0) = sum(allSign{j0}(:)  < 0)/length(allMat{j0}(:));
        stdRatio(j0+jstart,i0) = std(sum(allSign{j0} < 0,1)/length(allMat{j0}(:,1)));
        meanFmin(j0+jstart,i0) = - mean(allFmin{j0}(allFmin{j0} < 2));
        stdFmin(j0+jstart,i0) = std(allFmin{j0}(allFmin{j0} < 2));
    end       
end

% ==================================================================
% plot how inhibitory ratio changes with basal activity
% ==================================================================
figure
hold on
for i0 = 1:length(files)
    errorbar(R0',meanRatio(:,i0),stdRatio(:,i0),'o-','MarkerSize',...
    12,'Color',Bu(3 + 2*i0,:),'LineWidth',2,'CapSize',0)
end
plot([0;1],[0;1],'k--','LineWidth',2)
hold off
lg = legend(['N=',num2str(nOdor(1))],['N=',num2str(nOdor(2))],['N=',num2str(nOdor(3))],...
    ['N=',num2str(nOdor(4))]);
set(lg,'FontSize',20)
legend boxoff
box on
xlabel('$r_0$','Interpreter','latex','FontSize',28)
ylabel('inhibitory ratio','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_bothExciInhi_ratio_basal_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% ==================================================================
% plot how target function changes with basal activity
% ==================================================================
figure
hold on
for i0 = 1:length(files)
    errorbar(R0',meanFmin(:,i0),stdFmin(:,i0),'o-','MarkerSize',...
    12,'Color',Bu(3 + 2*i0,:),'LineWidth',2,'CapSize',0)
end
% plot([0;1],[0;1],'k--','LineWidth',2)
hold off
lg = legend(['N=',num2str(nOdor(1))],['N=',num2str(nOdor(2))],['N=',num2str(nOdor(3))],...
    ['N=',num2str(nOdor(4))]);
set(lg,'FontSize',20)
legend boxoff
box on
xlabel('$r_0$','Interpreter','latex','FontSize',28)
ylabel('differential entropy','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_bothExciInhi_fmin_basal_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% when paritial of the receptors have inhibitory interactions
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

%dfile = 'gcmi_partInhi_summData_N50sp2sig2_22-Jul-2018.mat';
dfile = 'gcmi_partInhi_summData_N50M20sp3sig2_07-Aug-2018.mat';

load(fullfile(dFolder,dfile))

% essential parameters, this can be get from the file nae
nOdor = 50;
nRecp = 20;
sig = 2;
sp = 2;
alp = 0.2;   % basal activity

% =============== plot example W =======================
inx = 1;
wi = reshape(allWi{2,1}(:,inx),[],nOdor);
we = reshape(allWe{2,1}(:,inx),[],nOdor);
w = [wi;we];
sign = reshape(Sign{1,1}(:,inx),[],nOdor);
%plot the heat map
figure
imagesc(log(w),[-8,6])
xlabel('odorant')
ylabel('receptor')
colormap(jet);
c = colorbar; 
c.FontSize = 16;

figNamePref = ['gcmi_partialInhi_N',num2str(nOdor),'M',num2str(nRecp),'sp',...
    num2str(sp),'sig',num2str(sig),'_hitmap_',date];
saveas(gcf,[sFolder,filesep,figNamePref,'.fig'])
print('-depsc',[sFolder,filesep,figNamePref,'.eps'])

% pd1 = fitdist(log(wi(sign <0)),'normal');
% pd2 = fitdist(log(wi(sign >0)),'normal');
% pd3 = fitdist(log(wi(sign <0)),'normal');
% 

% plot the histogram with Gaussian fit
[ha, pos] = tight_subplot(2,1,[0 0],[.15 .01],[.15 .03]);


% ======================================================
% plot a fitted histogram
% ======================================================
[bc1, edg1] = histcounts(log(wi(sign <0)),20);
[bc2, edg2] = histcounts(log(wi(sign >0)),20);

figure
hold on
Y1 = [bc1;bc1]/sum(bc1)/(edg1(2)- edg1(1)); %normalized, probability
ar1 = area(sort(edg1([1:end-1 2:end])),Y1(1:end));
ar1.FaceAlpha = 0.5;
ar1.FaceColor = Rd(4,:);
ar1.EdgeColor = Rd(8,:);
ar1.LineWidth = 2;


Y2 = [bc2;bc2]/sum(bc2)/(edg2(2)- edg2(1)); %normalized, probability
ar1 = area(sort(edg2([1:end-1 2:end])),Y2(1:end));
ar1.FaceAlpha = 0.5;
ar1.FaceColor = Bu(4,:);
ar1.EdgeColor = Bu(8,:);
ar1.LineWidth = 2;

hold off


% h1 = histfit(log(wi(sign <0)));
% h2= histfit(log(wi(sign >0)));
% % h3 = histfit(log(we(:)));
% hold off
% h1(1).Visible = 'off';h2(1).Visible = 'off';
% % h3(1).Visible = 'off';
% % h1(1).FaceColor = Bu(5,:);h2(1).FaceColor = Bu(5,:);h3(1).FaceColor = Bu(5,:);
% % h1(1).FaceAlpha = 0.4; h2(1).FaceAlpha = 0.4;h3(1).FaceAlpha = 0.4;
% 
% h1(2).Color = Bu(9,:);h2(2).Color = 'r';
% % h3(2).Color = 'k';
% h1(2).LineWidth = 3;h2(2).LineWidth = 3;
% h3(2).LineWidth = 3;

% legend(h1(2),'wi',h2(2),'wie',h3(2),'we')
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('frequency')
box on
figNamePref = ['gcmi_partialInhi_N',num2str(nOdor),'M',num2str(nRecp),'sp',num2str(sp),...
    'sig',num2str(sig),'_histW_',date];
saveas(gcf,[sFolder,filesep,figNamePref,'.fig'])
print('-depsc',[sFolder,filesep,figNamePref,'.eps'])


% ======================================================
% plot histogram of excitatory part
% ======================================================
figure
[bc1, edg1] = histcounts(log(we(:)),30);
Y1 = [bc1;bc1]/sum(bc1)/(edg1(2)- edg1(1)); %normalized, probability
ar1 = area(sort(edg1([1:end-1 2:end])),Y1(1:end));
ar1.FaceAlpha = 0.5;
ar1.FaceColor = Bu(4,:);
ar1.EdgeColor = Bu(8,:);
ar1.LineWidth = 2;


% histogram(log(we(:)),30)
xlabel('$\ln(w_e)$','interpreter','latex')
ylabel('frequency')
figNamePref = ['gcmi_partialInhi_N',num2str(nOdor),'M',num2str(nRecp),'sp',num2str(sp),...
    'sig',num2str(sig),'_histWe_',date];
saveas(gcf,[sFolder,filesep,figNamePref,'.fig'])
print('-depsc',[sFolder,filesep,figNamePref,'.eps'])

%% Analyze the new simulation data from Qianyi's simulation with both inihibtion and excitation
dFolder = '/Users/shan/Dropbox/olfactionProject/data/GcmiInhibit-06-13';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
allFile = dir([dFolder,filesep,'*.mat']);
allfiles = {allFile.name}';

nOdor = [20,30,40,50];     %need to specify this part for different data set
str_marker = '2018-06-22';   
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(allfiles,str,'once'));
files = allfiles(FIND(str_marker));

% default graphics settings
defaultGraphicsSetttings

%initilze the variables
% allRatio = cell(length(nRecp),1);
% stdRatio = cell(length(nRecp),1);
% 
% meanFmin = cell(length(nRecp),1);
% stdFmin = cell(length(nRecp),1);

allRatio = zeros(length(files),1);
stdRatio = zeros(length(files),1);
meanFmin = zeros(length(files),1);
stdFmin = zeros(length(files),1);
allr0 = zeros(length(files),1);
allN = zeros(length(files),1);
% exciEle = zeros(length(files),2); %mean and std
% inhiEle = zeros(length(files),2); %mean and std


for i0 = 1:length(files)
    s1 = '(?<= *N)[\d.]+(?=_)';
    s2 = '(?<= *basal)[\d.]+(?=_)';
    allN(i0) = str2num(char(regexp(files{i0},s1,'match')));
    allr0(i0) = str2num(char(regexp(files{i0},s2,'match')));
    
    load(fullfile(dFolder,files{i0}))
    allRatio(i0) = sum(allSign(:) < 0)/length(allMat(:));
    stdRatio(i0) = std(sum(allSign < 0,1)/length(allMat(:,1)));
    meanFmin(i0) = - mean(allfmin(allfmin < 2));  % here 2 is used to get rid of outliner
    stdFmin(i0) = std(allfmin(allfmin < 2));
%     exciEle(i0,:) =  [mean(log(allMat(allSign > 0))),std(log(allMat(allSign > 0)))];
%     inhiEle(i0,:) =  [mean(log(allMat(allSign < 0))),std(log(allMat(allSign < 0)))];
end

% [uniq_r0,rInx] = uniq(allr0);
uniq_r0 = [0.05,0.1,0.15,0.2:0.1:0.8,0.85,0.9,0.95];
NR = length(uniq_r0);
NO = length(nOdor);
summaryData  = struct('allRatio',nan(NR,NO),'stdRatio',nan(NR,NO),'meanFmin',...
    nan(NR,NO),'stdFmin',nan(NR,NO),'nRecp',nOdor,'allR0',uniq_r0);
for i0 = 1:length(nOdor)
    inx =find(allN == nOdor(i0));
    r0 = allr0(inx);
%     [~,order] = sort(round(r0*100));
    [~,r0_inx] = ismember(round(r0*100),round(uniq_r0*100));
    summaryData.allRatio(r0_inx,i0) = allRatio(inx);
    summaryData.stdRatio(r0_inx,i0) = stdRatio(inx);
    summaryData.meanFmin(r0_inx,i0) = meanFmin(inx);
    summaryData.stdFmin(r0_inx,i0) = stdFmin(inx);
end

%save the data
sfName = ['gcmi_inhiActi_diffBasal_diffN_QY_',date,'.mat'];
save(fullfile(outFolder,sfName),'-struct','summaryData')

% ===============================================================
% plot the how differential entropy change with basal activity
% ===============================================================
figure
hold on
for i0 = 1:length(nOdor)
    errorbar(uniq_r0',summaryData.meanFmin(:,i0),summaryData.stdFmin(:,i0),...
        'o-','MarkerSize',12,'Color',Bu(3 + 2*i0,:),'LineWidth',2,'CapSize',0)
end
hold off
lg = legend(['N=',num2str(nOdor(1))],['N=',num2str(nOdor(2))],['N=',num2str(nOdor(3))],...
    ['N=',num2str(nOdor(4))],'Location','northwest');
set(lg,'FontSize',20)
legend boxoff
box on
xlabel('$r_0$','Interpreter','latex','FontSize',28)
ylabel('differential entropy','FontSize',28)
set(gca,'XTick',0:0.2:1,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_bothExciInhi_fmin_basal_QY',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ===============================================================
% plot the how inhibitory ratio change with basal activity
% ===============================================================
figure
hold on
for i0 = 1:length(nOdor)
    errorbar(uniq_r0',summaryData.allRatio(:,i0),summaryData.stdRatio(:,i0),...
        'o-','MarkerSize',12,'Color',Bu(3 + 2*i0,:),'LineWidth',2,'CapSize',0)
end
plot([0;1],[0;1],'k--','LineWidth',2)
hold off
lg = legend(['N=',num2str(nOdor(1))],['N=',num2str(nOdor(2))],['N=',num2str(nOdor(3))],...
    ['N=',num2str(nOdor(4))],'Location','northwest');
set(lg,'FontSize',20)
legend boxoff
box on
xlabel('$r_0$','Interpreter','latex','FontSize',28)
ylabel('inhibitory ratio','FontSize',28)
set(gca,'XTick',0:0.2:1,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_bothExciInhi_ratio_basal_QY',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% Analyze my new Gcmi inhibition data 10/13/2018
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_inhi/N50M10S2sig2_1013';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir([dFolder,filesep,'*.mat']);
files = {allFile.name}';

N = 50;
M = 10;
sig = 10;
sp = 2;

R0 = [0.05,0.08,0.1,0.12,0.14,0.16,0.18,0.25,0.3,0.4,0.5,0.6,0.7,0.82,0.84,0.88,0.9,0.95];  % the basal activity

meanRatio = nan(length(R0),1);
stdRatio = nan(length(R0),1);
meanFmin = nan(length(R0),1);
stdFmin = nan(length(R0),1);
for i0 = 1:length(files)
    s2 = '(?<= *_alp)[\d.]+(?=_)';
    r0 = str2num(char(regexp(files{i0},s2,'match')));
    inx = find(100*R0 == 100*r0);

    load(fullfile(dFolder,files{i0}));
    
    meanRatio(inx) = sum(allSign(:) < 0)/length(allMat(:));
    stdRatio(inx) = std(sum(allSign < 0,1)/length(allMat(:,1)));
    meanFmin(inx) = - mean(allfmin(allfmin < 2));  % here 2 is used to get rid of outliner
    stdFmin(inx) = std(allfmin(allfmin < 2));    
end
%save the data
sfName = ['gcmi_inhi_N50M10sp2_Acti_diffBasal__',date,'.mat'];
save(fullfile(outFolder,sfName),'meanRatio','stdRatio','meanFmin','stdFmin',...
    'R0','N','M','sig','sp')

% ==========================================
% plot how inhibitory fraction change with r0
% ==========================================
figure
hold on
errorbar(R0',meanRatio,stdRatio,'o-','MarkerSize',12,'MarkerFaceColor',Bu(9,:),...
    'MarkerEdgeColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
plot([0;1],[0;1],'k--','LineWidth',2)
hold off
box on
xlabel('$r_0$','Interpreter','latex','FontSize',28)
ylabel('inhibitory fraction','FontSize',28)
set(gca,'XTick',0:0.2:1,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_bothExciInhi_ratio_basal_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% ==========================================
% plot how fmin change with r0
% ==========================================
figure
hold on
errorbar(R0',meanFmin,stdFmin,'o-','MarkerSize',12,'MarkerFaceColor',Bu(9,:),...
    'MarkerEdgeColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
hold off
box on
xlabel('$r_0$','Interpreter','latex','FontSize',28)
ylabel('differential entropy','FontSize',28)
set(gca,'XTick',0:0.2:1,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_bothExciInhi_fmin_basal_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


%% N = 1, M --> infinity, comparison with theory

close all

% data folder and output folder
dFolder = '/Users/shan/Dropbox/olfactionProject/data/oneByM01142018';
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
        allWmean(j0,i0) = mean(temp(:));
        allWstd(j0,i0) = std(temp(abs(temp) <=4));
    end
end

% theoretical std of wij
h = 1;
[orderSig,inx] = sort(allSig);
theoStd = sqrt(max(orderSig.^2 - 1/h^2,0));  % when sigma_c is smaller than 1

%plot the and compare with theory
% myColor = brewermap(11,'Blues');

% ===== differential entropy =============
% ========================================
figure
hold on
largeMeanStd = zeros(2,NUM); %store the mean and std of last data, very large M
for i0 = 1:NUM
%     errorbar(allM',allFmin(:,inx(i0)),allFminStd(:,inx(i0)),'Color',Gr(2+i0,:),...
%         'LineWidth',3)
    largeMeanStd(1,inx(i0)) = allFmin(end,inx(i0));
    largeMeanStd(2,inx(i0)) = allFminStd(end,inx(i0));
end
hold off
xlabel('M','FontSize',28)
ylabel('differential entropy','FontSize',28)
set(gca,'XTick',[1,10,100],'XScale','log','FontSize',24,'LineWidth',1.5)
figNamePref = 'oneByMdiffEntropySimu';
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ==========plot differential entropy as a function of sigma_c ==========
% ======================================================================
figure 
errorbar(orderSig,largeMeanStd(1,inx),largeMeanStd(2,inx),'o-','MarkerSize',...
    12,'MarkerFaceColor',Bu(8,:),'Color',Bu(8,:),'LineWidth',3)
legend('M=200','Location','northwest')
legend boxoff

xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('differential entropy','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','log')
figNamePref = 'oneByMdiffEntropyLargeM';
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ========= compare the wij with theory =====================
% ===========================================================
figure
hold on
allStdOrder = sort(allWstd(25,:).*allSig);
plot(theoStd,allStdOrder,'o','MarkerEdgeColor',Bu(9,:),...
    'MarkerFaceColor',Bu(9,:),'LineWidth',2,'MarkerSize',12)
plot([0,10],[0,10],'k--','LineWidth',2)
hold off
lg = legend('simulation','theory','Location','northwest');
legend boxoff

xlim([0 11])
ylim([0 11])
xlabel('$\sqrt{\sigma_c^2 - 1/h^2}$','Interpreter','latex','FontSize',28)
ylabel('$\sigma_w$','Interpreter','latex','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5)
figNamePref = 'oneByMCompSimuTheory';
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


S = 0:0.1:10;
theorSig = sqrt(max(S.^2 - 1/h^2,0));
figure
hold on
plot(S,theorSig,'k-','LineWidth',2)
plot(sort(allSig),allStdOrder,'o','MarkerEdgeColor',Bu(9,:),...
    'MarkerFaceColor',Bu(9,:),'MarkerSize',12)
lg = legend('theory','simulation','Location','northwest');
legend boxoff
xlim([-0.2 11])
ylim([-0.2,11])
xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('$\sigma(W)$','Interpreter','latex','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1)
figNamePref = 'oneByMsigWSimu';
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

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
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

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


% =====================================================================
% the average position interval, to see if they are uniform in the large M
% limit
% ======================================================================
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
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


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

figFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
fileNameEps = [figFolder,filesep,'normpdf_cdf.eps'];
fileNameFig = [figFolder,filesep,'normpdf_cdf.fig'];
print('-depsc',fileNameEps)
saveas(gcf,fileNameFig)

%% N = 2, M --> infinity, comparison with theory

% for h = 1
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
fName = 'N2diffSigSummary.mat';
load(fullfile(dFolder,fName))

%% N = 2, M --> infinity, how sparsity of w changes with input sigma
% plot the differential entropy as a function of sigma_c
% dFolder = '/Users/shan/Dropbox/olfactionProject/data/N2difSig/DifSigData';
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

file = 'summData_0418_Regul.mat';
load(fullfile(dFolder,file))
h = 1;  %Hill Coefficient
% allFile = dir([dFolder,filesep,'*.mat']);
% allfiles = {allFile.name}';
% allSig = zeros(length(allfiles),1);
% spStat = zeros(length(allfiles),2);  % mean and std of sparisty
% 
% for i0 = 1:length(allfiles)
%     str = '(?<= *sig)[\d.]+(?=_)';
%     allSig(i0) = str2num(char(regexp(allfiles{i0},str,'match')));
%     load(fullfile(dFolder,allfiles{i0}))
%     
%     temp = allMat{end};
%     if size(temp,1) ~= 400
%         error('the odorant is not 200')
%     else
%         spStat(i0,1) = sum(temp(:) > -3)/length(temp(:));
%         spStat(i0,2) = std(sum(temp > -3,1)/size(temp,1));
%         
%     end
% end
% ========================================================================
% plot the how sparsity change with sigma_c
% ========================================================================
inx = [1,2,4];
figure
hold on
for i0 = 1:length(inx)
errorbar(allSig',spMean(inx(i0),:),spStd(inx(i0),:),'o-','MarkerSize',...
    12,'MarkerFaceColor',Bu(2+3*i0,:),'Color',Bu(2+3*i0,:),'CapSize',0,'LineWidth',2)
end
hold off
% legend('M=200','Location','northwest')
legend('M=60','M=100','M=200','Location','northwest')
legend boxoff

xlabel('$\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('sparsity of W','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = 'N2M_spW_sigc_h1';
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ========================================================================
% plot the how fmin depends on input sigma
% ========================================================================
figure
hold on
for i0 = 1:length(inx)
errorbar(allSig',allFmean(inx(i0),:),allFstd(inx(i0),:),'o-','MarkerSize',...
    12,'MarkerFaceColor',Bu(2+3*i0,:),'Color',Bu(2+3*i0,:),'CapSize',0,'LineWidth',2)
end
hold off
% legend('M=200','Location','northwest')
legend('M=60','M=100','M=200','Location','northwest')
legend boxoff

xlabel('$\sigma_c$','Interpreter','latex')
ylabel('differential entropy')
set(gca,'XScale','linear')
figNamePref = 'N2M_fmin_sigc_h1';
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% N -->, M--> infity, but n is finite, a mean field theory

%% optimal parameterized distribution
dFolder = '/Users/shan/Dropbox/olfactionProject/data/intEqnSolve_parameterized';
sigall = [1.0:0.1:3.1,3.3:0.1:4.4,4.6:0.1:4.9];
for i0 = 1:length(sigall)
  fName = strcat('Int_sig',num2str(sigall(i0)),'.mat');
  load(fullfile(dFolder,fName));  
  meanvalue(i0) = x(1);
  stdvalue(i0) = x(2);
  sparsity(i0) = x(3);
end
subplot(1,3,1)
plot(sigall,sparsity,'linewidth',3);
xlabel('\sigma');
ylabel('ln(w) sparsity');
set(gca,'linewidth',1.5,'fontsize',20)
subplot(1,3,2)
plot(sigall,meanvalue,'linewidth',3);
xlabel('\sigma');
ylabel('ln(w) \mu');
set(gca,'linewidth',1.5,'fontsize',20);
subplot(1,3,3)
plot(sigall,stdvalue,'linewidth',3);
xlabel('\sigma');
ylabel('ln(w) \sigma');
set(gca,'linewidth',1.5,'fontsize',20)

%% Reconstruction or decoding
% data from Qianyi's simulation
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding/recons_H100_N20_R9_S2_sig2_noise0_2018-09-27';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dFolder,filesep,'*.mat'));
files = {allFile.name}';

% set the parameters
N = 20;         % number of odorants
M = 9;         % number of receptors
sp = 2;         % sparsity of input
sig = 2;        % std of input
noiseSig =0.0001; % noise level 
H = 50;        %hidden layer size
% allSpInput = 0.05:0.05:0.95;
allSpInput = 0.1:0.1:1;

allLayer = 1:1:3;  %number of hidden layers;
% meanTrainPf = zeros(length(allSp),length(allLayer));  %mean and std of best trianing performance
% stdTrainPf = zeros(length(allSp),length(allLayer));  %mean and std of best trianing performance

% meanTestPf = zeros(length(allSp),length(allLayer));  %mean and std of best trianing performance
% stdTestPf = zeros(length(allSp),length(allLayer));   %mean and std of best trianing performance

allTrain = zeros(length(files),length(allSpInput),length(allLayer));
allTest = zeros(length(files),length(allSpInput),length(allLayer));

for i0 = 1:length(files)
    load(fullfile(dFolder,files{i0}))
    for j0 = 1:length(allSpInput)
        for k0 = 1:length(allLayer)
            allTrain(i0,j0,k0) = tr_all{j0,k0}.best_perf;
            allTest(i0,j0,k0) = tr_all{j0,k0}.best_tperf;
        end
    end
end

meanTrainPf = squeeze(mean(allTrain,1));
stdTrainPf = squeeze(std(allTrain,0,1));

meanTestPf = squeeze(mean(allTest,1));
stdTestPf = squeeze(std(allTest,0,1));

% log scale
logMeanTestPf = squeeze(mean(log10(allTest),1));
logStdTestPf = squeeze(std(log10(allTest),0,1));

minLogPf = min(logMeanTestPf,[],1);
% ===============================================================
% plot the sparsity-dependent reconstruction error
% ==============================================================
figure
hold on
for i0 = 1:length(allLayer) 
    errorbar(allSpInput',meanTestPf(:,i0),stdTestPf(:,i0),'o-','MarkerSize',...
        12,'MarkerFaceColor',Bu(2+i0*3,:),'Color',Bu(2+i0*3,:),'LineWidth',2,...
        'CapSize',0)
end
hold off
lg = legend([num2str(allLayer(1)),' hidden layer'],[num2str(allLayer(2)),' hidden layer'],...
    [num2str(allLayer(3)),' hidden layer']);
set(lg,'FontSize',20)
legend boxoff
ylim([5e-2,0.15])
xlabel('sparsity of W')
ylabel('MSE')
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear','Yscale','linear')
figNamePref = ['ReconstructMSE_N',num2str(N),'_M',num2str(M),'_sp',num2str(sp),...
    '_noise',num2str(noiseSig),'_diffSpW_linear',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% log scale
sInx = [1,3]; % only select two 
figure
hold on
for i0 = 1:length(sInx) 
    errorbar(allSpInput',logMeanTestPf(:,sInx(i0)),logStdTestPf(:,sInx(i0)),'o-','MarkerSize',...
        12,'MarkerFaceColor',Bu(3+sInx(i0)*2,:),'Color',Bu(3+sInx(i0)*2,:),'LineWidth',2,...
        'CapSize',0,'LineWidth',2)

end
hold off
% lg = legend([num2str(allLayer(1)),' hidden layer'],[num2str(allLayer(2)),' hidden layer'],...
%     [num2str(allLayer(3)),' hidden layer'],[num2str(allLayer(4)),' hidden layer']);
lg = legend([num2str(allLayer(sInx(1))),' hidden layer'],[num2str(allLayer(sInx(2))),' hidden layer']);
set(lg,'FontSize',20)
legend boxoff
box on

xlabel('sparsity of W')
ylabel('$\log10$(MSE)','Interpreter','latex')
set(gca,'FontSize',24,'LineWidth',1.5)
figNamePref = ['ReconstructMSE_N20M9_sp2_sig2_diffSpW_log',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% only plot one hidden layer
figure
errorbar(allSpInput',logMeanTestPf(:,1),logStdTestPf(:,1),'o-','MarkerSize',...
        12,'MarkerFaceColor',Bu(10,:),'Color',Bu(10,:),'LineWidth',2,...
    'CapSize',0)

lg = legend([num2str(allLayer(1)),' hidden layer'],[num2str(allLayer(2)),' hidden layer'],...
    [num2str(allLayer(3)),' hidden layer'],[num2str(allLayer(4)),' hidden layer']);
set(lg,'FontSize',20)
legend boxoff

xlabel('sparsity of W')
ylabel('$\log10(MSE)$','Interpreter','latex')
set(gca,'FontSize',24,'LineWidth',1.5)
figNamePref = ['Reconstruct_MSE_N20M9_sp2_sig2_diffSpW_oneHidden_log',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% ======= How the best performance depends on number of hidden layers ====
figure
plot((1:1:4),minLogPf,'o-','MarkerSize',12,'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:))
ylim([-3,-2.5]);
xlabel('number of hidden layer')
ylabel('$\log(\rm{MSE})$','interpreter','latex')
figNamePref = ['Recon_MSE_N20M9_sp2_sig2_diffSpW_diffLayer_log',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% Reconstruction with new data, 5 different layers
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
dName = 'recons_N50M10sp2_noise0.01H100thd0.13534_24-Sep-2018.mat';
load(fullfile(dFolder,dName))

% basic parameters
N = 50;
M = 10;
sp = 2;
sig = 2;
noiseSig = 0.01;
H = 100;
thd = exp(-sig);  %threshold to be detected, default 1

% ==================================================
% plot the sparsity with respect to r0
% ==================================================
meanMSE = mean(allMSE,3);
stdMSE = std(allMSE,0,3);
minMSE = min(meanMSE,[],1);

figure
hold on
% errorbar(allSp'*ones(1,5),meanMSE,stdMSE)
for i0 = 1:size(meanMSE,2)
    plot(allSp',meanMSE(:,i0),'Color',Bu(1+2*i0,:))
end
hold off
xlabel('$r_0$','interpreter','latex')
ylabel('$MSE$','interpreter','latex')
figNamePref = ['Recon_MSE_N',num2str(N),'M',num2str(M),'_sp',num2str(sp),...
    '_noiseSig',num2str(noiseSig),'_thd',nu2str(thd),'_H',num2str(H),'_diffSpW_',date];
% figNamePref = ['Recon_MSE_N20M10_sp2_sig2_diffSpW_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])



% ==================================================
% plot the optimal performance vs number of hidden layers
% ==================================================
figure
plot(1:1:5,minMSE,'o-','MarkerSize',12,'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:))
ylim([0.03,0.05])
xlabel('number of hidden layer')
ylabel('$\rm{MSE}$','interpreter','latex')
% figNamePref = ['Recon_MSE_N20M10_sp2_sig2_diffSpW_hiddenLayers',date];
figNamePref = ['Recon_MSE_N',num2str(N),'M',num2str(M),'_sp',num2str(sp),...
    '_noiseSig',num2str(noiseSig),'_thd',nu2str(thd),'_H',num2str(H),'_hiddenLayers_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% Reconstruction using tensorflow, data from Qianyi, how performance change with sparsity
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
dFolder = '/Users/shan/Dropbox/olfactionProject/data/reconstruction_nn';
% outFolder = '/Users/shan/Dropbox/olfactionProject/data/reconstruction_nn';

dName1 = 'results_difsp_100*20.mat';

load(fullfile(dFolder,dName1))
N = 100;
M = 20;
sp = 2;
sig = 2;
noisSig = 0.05;

allSp = 1:-0.1:0.1;

randRef = loss_fcn(end,:);
meanLoss = mean(loss_fcn(1:end-1,:),2);
stdLoss = std(loss_fcn(1:end-1,:),0,2);

logMeanLoss = mean(log10(loss_fcn(1:end-1,:)),2);
logStdLoss = std(log10(loss_fcn(1:end-1,:)),0,2);

% ==============================
% loss function at linear scale
% ==============================
figure
hold on
errorbar(allSp',meanLoss,stdLoss,'o-','MarkerSize',10,'MarkerFaceColor',...
    Bu(9,:),'Color',Bu(9,:),'LineWidth',1.5,'CapSize',0)
errorbar(0,mean(randRef),std(randRef),'o','MarkerSize',10,'MarkerFaceColor',...
    RdBu(3,:),'Color',RdBu(3,:),'LineWidth',1.5,'CapSize',0)
hold off
% ylim([0.7,1])
set(gca,'Xtick',0:0.5:1)
xlabel('$\rho_w$','interpreter','latex')
ylabel('loss function')
figNamePref = ['recons_loss_tensorflow_N',num2str(N),'M',num2str(M),...
    'sp',num2str(sp),'_nSig',num2str(noisSig),'_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% ==============================
% loss function at log scale
% ==============================
figure
hold on
errorbar(allSp',logMeanLoss,logStdLoss,'o-','MarkerSize',10,'MarkerFaceColor',...
    Bu(9,:),'Color',Bu(9,:),'LineWidth',1.5,'CapSize',0)
% errorbar(0,mean(log10(randRef)),std(log10(randRef)),'o','MarkerSize',10,'MarkerFaceColor',...
%     RdBu(3,:),'Color',RdBu(3,:),'LineWidth',1.5,'CapSize',0)
hold off
lg = legend(['N=',num2str(N), ',M =',num2str(M), ',n =',num2str(sp),',\sigma_c = ',num2str(sig)]);
set(lg,'FontSize',16)
legend boxoff
box on
% ylim([0.7,1])
set(gca,'Xtick',0:0.5:1)
xlabel('$\rho_w$','interpreter','latex')
ylabel('$\log10$(loss function)','interpreter','latex')
figNamePref = ['recons_loss_tensorflow_N',num2str(N),'M',num2str(M),...
    'sp',num2str(sp),'_nSig',num2str(noisSig),'_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ========================================================
% Scatter plot to show the performance of reonstruction
% =========================================================

% first select the indices
ix = 10; iy = 3;
x0 = squeeze(original(ix,iy,:,:));
x_hat = squeeze(prediction(ix,iy,:,:));

% histogram and scatter plot, using tight plot
figureSize = [0 0 12 4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(1,3,[.05 .08],[.25 .05],[.1 .03]);

axes(ha(1))
temp = x0; temp(temp ==0) = -100;
imagesc(temp,[-6,6])
xlabel('odorant index')
ylabel('trial')

axes(ha(2))
temp = x_hat;temp(temp == 0) = -100;
imagesc(temp,[-6,6])
xlabel('odorant index')
ylabel('trial')

%scatter plot of the 
axes(ha(3));
nonZeorX0 = zeros(size(x0,1),sp);
nonZeorXh = zeros(size(x0,1),sp);
for i0 = 1:size(x0,1)
    nonZeorX0(i0,:) = sort(x0(i0,x0(i0,:) ~=0));
    nonZeorXh(i0,:) = sort(x_hat(i0,x_hat(i0,:) ~=0));
end
scatter(nonZeorX0(:,1),nonZeorXh(:,1),20,Bu(9,:),'filled')
hold on
scatter(nonZeorX0(:,2),nonZeorXh(:,2),20,Bu(9,:),'filled')
xl = ha(3).XLim;
plot(xl',xl','--','LineWidth',2,'Color',Gr(8,:))
hold off
xlabel('$x$','interpreter','latex')
ylabel('$\hat{x}$','interpreter','latex')

figNamePref = ['recons_tensorflow_N',num2str(N),'M',num2str(M),...
    'sp',num2str(sp),'_nSig',num2str(noisSig),'_sigw_',num2str(allSp(ix)),'_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% axes(ha(4))
% scatter(nonZeorX0(:,2),nonZeorXh(:,2),20,Bu(9,:),'filled')


% =======================================================
% scatter plot to compare the reconstruction performance
% ==============================================
% dFolder = '/Users/shan/Documents/machine learning/olfactionEnco';
dFolder = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/code/decoding/data/reconsFig6';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% dNames = {'loss_RandRhoW_N100M20noise0.05_sp0.9.mat',...
%             'loss_RandRhoW_N100M20noise0.05_sp0.45.mat',...
%             'loss_RandRhoW_N100M20noise0.05_sp0.05.mat'};
        
dNames = {'loss_hardWire_100M20noise0.05_sp0.9_L2.mat',...
            'loss_hardWire_100M20noise0.05_sp0.4_L2.mat',...
            'loss_hardWire_100M20noise0.05_sp0.05_L2.mat'};
        
% dNames = {'loss_hardWire_100M20noise0.1_sp0.9000000000000001_L2.mat',...
%             'loss_hardWire_100M20noise0.1_sp0.4_L2.mat',...
%             'loss_hardWire_100M20noise0.1_sp0.05_L2.mat'};
        
% parameters
N = 100; M = 20; sp = 2; sig = 2; 
noise = 0.05;
allSp = [0.1, 0.6, 0.95];

% set graph size
figureSize = [0 0 13 10];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% tight plot
[ha, pos] = tight_subplot(3,3,[.05 .1],[.1 .02],[.1 .03]);

inx = 501:600;  %select partial data

for i0 = 1:3
    load(fullfile(dFolder,dNames{i0}))
    
    axes(ha(3*i0-2))
    temp1 = original(inx,:); temp1(temp1 ==0) = -100;
    imagesc(temp1,[-6,6])
    
    ylabel('trial')
    if i0==3 
        xlabel('odorant index')
    else
        set(gca,'XTick',[])
    end

    axes(ha(3*i0-1))
    temp2 = prediction(inx,:);temp2(temp2 == 0) = -100;
    imagesc(temp2,[-6,6])
%     xlabel('odorant index')
    ylabel('trial')
    if i0==3
        xlabel('odorant index')
    else
        set(gca,'XTick',[])
    end   
    
    axes(ha(3*i0));
    temp1 = original;
    temp2 = prediction;
    nonZeorX0 = zeros(size(temp1,1),sp);
    nonZeorXh = zeros(size(temp2,1),sp);
    for j0 = 1:size(temp1,1)
        nonZeorX0(j0,:) = sort(temp1(j0,temp1(j0,:) ~=0));
        nonZeorXh(j0,:) = sort(temp2(j0,temp2(j0,:) ~=0));
    end
    scatter(nonZeorX0(:,1),nonZeorXh(:,1),10,Bu(8,:),'filled')
    hold on
    scatter(nonZeorX0(:,2),nonZeorXh(:,2),10,Bu(8,:),'filled')
    xlim([-5.5,7])
    ylim([-5.5,7])
    xl = ha(3*i0).XLim;
    plot(xl',xl','--','LineWidth',2,'Color',Gr(8,:))
    hold off
    box on
    ylabel('$\ln(\hat{c})$','interpreter','latex')
    if i0==3
        xlabel('$\ln(c)$','interpreter','latex')
    else
        set(gca,'XTick',[])
    end
end

figNamePref = ['recons_tensorflow_N',num2str(N),'M',num2str(M),...
    'sp',num2str(sp),'_nSig',num2str(noise),'_scatter_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% ===========
% plot a scatter in a row
% ===========
figureSize = [0 0 12 4.3];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% tight plot
[ha, pos] = tight_subplot(1,3,[.05 .03],[.25 .08],[.1 .03]);
for i0 = 1:3
    load(fullfile(dFolder,dNames{i0}))
    
    axes(ha(i0));
    temp1 = original;
    temp2 = prediction;
    nonZeorX0 = zeros(size(temp1,1),sp);
    nonZeorXh = zeros(size(temp2,1),sp);
    for j0 = 1:size(temp1,1)
        nonZeorX0(j0,:) = sort(temp1(j0,temp1(j0,:) ~=0));
        nonZeorXh(j0,:) = sort(temp2(j0,temp2(j0,:) ~=0));
    end
    scatter(nonZeorX0(:,1),nonZeorXh(:,1),10,Bu(7,:),'filled')
    hold on
    scatter(nonZeorX0(:,2),nonZeorXh(:,2),10,Bu(7,:),'filled')
    lg = legend(['Error = ',num2str(round(testLost,2))],'Location','northwest');
    set(lg,'FontSize',16)
    legend boxoff
    xlim([-5.5,7])
    ylim([-5.5,7])
    xl = ha(i0).XLim;
    plot(xl',xl','--','LineWidth',1.5,'Color',Gr(8,:))
    hold off
    box on
    xlabel('$\ln(c)$','interpreter','latex')
    if i0==1
        ylabel('$\ln(\hat{c})$','interpreter','latex')
    else
        set(gca,'YTick',[])
    end
end
figNamePref = ['recons_tensorflow_N',num2str(N),'M',num2str(M),...
    'sp',num2str(sp),'_nSig',num2str(noise),'_scatter_exp_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])



% ====================================
% only compare reconstruction heatmap
% ====================================
% set graph size
figureSize = [0 0 13 6];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% tight plot
[ha, pos] = tight_subplot(2,3,[.05 .08],[.15 .02],[.1 .08]);

inx = 501:600;  %select partial data

for i0 = 1:3
    load(fullfile(dFolder,dNames{i0}))
    
    axes(ha(i0))
    temp1 = original(inx,:); temp1(temp1 ==0) = -100;
    imagesc(temp1,[-6,6])
    set(gca,'XTick',[])
    ylabel('trial')
%   

    axes(ha(i0+3))
    temp2 = prediction(inx,:);temp2(temp2 == 0) = -100;
    imagesc(temp2,[-6,6])
    ylabel('trial')
    xlabel('odorant index')
%     if i0==2
%         xlabel('odorant index')
%     else
%         set(gca,'XTick',[])
%     end   
    
%     axes(ha(3*i0));
%     temp1 = original;
%     temp2 = prediction;
%     nonZeorX0 = zeros(size(temp1,1),sp);
%     nonZeorXh = zeros(size(temp2,1),sp);
%     for j0 = 1:size(temp1,1)
%         nonZeorX0(j0,:) = sort(temp1(j0,temp1(j0,:) ~=0));
%         nonZeorXh(j0,:) = sort(temp2(j0,temp2(j0,:) ~=0));
%     end
%     scatter(nonZeorX0(:,1),nonZeorXh(:,1),10,Bu(8,:),'filled')
%     hold on
%     scatter(nonZeorX0(:,2),nonZeorXh(:,2),10,Bu(8,:),'filled')
%     xlim([-5.5,7])
%     ylim([-5.5,7])
%     xl = ha(3*i0).XLim;
%     plot(xl',xl','--','LineWidth',2,'Color',Gr(8,:))
%     hold off
%     box on
%     ylabel('$\ln(\hat{c})$','interpreter','latex')
%     if i0==3
%         xlabel('$\ln(c)$','interpreter','latex')
%     else
%         set(gca,'XTick',[])
%     end
end
figNamePref = ['recons_tensorflow_N',num2str(N),'M',num2str(M),...
    'sp',num2str(sp),'_nSig',num2str(noise),'_heatMap_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% Reconstruction within Inhibition, comparison

dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
fName = 'recons_Inhi_tensor_N50M10sp2sig2ns0.01_loss_24-Dec-2018New.mat';
refFile = 'recons_tensor_N50M10sp2sig2ns0.01_loss_24-Dec-2018.mat'; % reference
summData = load(fullfile(dFolder,fName));

% load the excitation-only data
load(fullfile(dFolder,refFile))
aveError = mean(allTestLoss,2);
stdError = std(allTestLoss,0,2);
aveMin = min(aveError);



nOdor = 50;
nRecp = 10;
noisSig = 0.05;
H = 200;

figure
hold on
errorbar(summData.allSp',nanmean(summData.allTestLoss,2),nanstd(summData.allTestLoss,0,2),'o-','MarkerSize',10,'MarkerFaceColor',...
    Bu(9,:),'Color',Bu(9,:),'LineWidth',1.5,'CapSize',0)
% set(lg,'FontSize',16)
plot([0;1],[aveMin;aveMin],'k--')
ylim([0.3,0.55])
box on
set(gca,'Xtick',0:0.2:1)
xlabel('$r_0$','interpreter','latex')
ylabel('reconstruction error')
figNamePref = ['class_Error_N',num2str(nOdor),'M',num2str(nRecp),...
    'sp',num2str(sp),'_nSig',num2str(noisSig),'_group',num2str(group),'_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% plot the error v.s input sparsity
figure
errorbar(1-allSp',aveError,stdError,'o-','MarkerSize',10,'MarkerFaceColor',...
    Bu(9,:),'Color',Bu(9,:),'LineWidth',1.5,'CapSize',0)
xlabel('$\rho_w$','interpreter','latex')
ylabel('reconstruction error')

%% classification
% this section plot how classification error changes with input 
% sparsity


dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

nOdor = 50;
nRecp = 10;
noisSig = 0.01;

NP = 1e3;   %number of different patterns

% allSpInput = [2,3,4,5];
allSpInput = 3;

allSpW = 0.05:0.05:1;

meanPr = zeros(length(allSpW),length(allSpInput));
stdPr = zeros(length(allSpW),length(allSpInput));

for i0 = 1:length(allSpInput)
%     fName = ['classify_N20M9sp',num2str(allSpInput(i0)),'_noiseSig0.01_28-Jun-2018.mat'];
    fName = ['classify_N50_R10_S',num2str(allSpInput(i0)),'_noiseSig0.01_nType3_2018-07-18'];
    load(fullfile(dFolder,fName))
    [Nsp,Rp] = size(allCrossEntr);  %number of sparsity and repeats
    meanPr(:,i0) = mean(exp(-allCrossEntr),2);
    stdPr(:,i0) = std(exp(-allCrossEntr),0,2);
end
% fName = 'classify_N20M9sp3_noiseSig0.01_28-Jun-2018.mat';
% fName = 'classify_N20M9sp4_noiseSig0.01_28-Jun-2018.mat';
% fName = 'classify_N20M9sp5_noiseSig0.01_28-Jun-2018.mat';


% allTestPr = zeros(Nsp,Rp);
%     for j0 = 1:Rp
%         temp = summData(:,j0).allTpr;
%         allTestPr(:,j0) = temp(:,1);
%     end
    
% meanPr = mean(log10(allCrossEntr),2);
% stdPr = std(log10(allCrossEntr),0,2);


% meanPr = mean(allTestPr,2);
% stdPr = std(allTestPr,0,2);

% plot the figure
inx = [2,3,4]; %only plot three
figure
hold on
for i0 = 1:length(inx)
errorbar(allSpW',meanPr(:,inx(i0)),stdPr(:,inx(i0)),'o-','MarkerSize',10,'MarkerFaceColor',...
    Bu(2 + 3*i0,:),'Color',Bu(2 + 3*i0,:),'LineWidth',1.5,'CapSize',0)
end
hold off
lg = legend(['sp = ',num2str(allSpInput(inx(1)))],['sp = ',num2str(allSpInput(inx(3)))],...
    ['sp = ',num2str(allSpInput(inx(3)))],'Location','southeast');

set(lg,'FontSize',16)
legend boxoff
box on
xlabel('sparsity of W')
% ylabel('$\log10$(cross entropy)','Interpreter','latex')
ylabel('accuracy')

set(gca,'FontSize',24,'LineWidth',1.5,'XTick',0:0.2:1)
figNamePref = ['classification_crossEntr_N',num2str(nOdor),'M',num2str(nRecp),...
    'diffSp','_nSig',num2str(noisSig),'_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])



% ================================================================
% Classification with different noise around centroid, here noise is 0.1
% ================================================================
file = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/decoding/classify_LDA_N50M12_HS500_sp2_group2_npattern100_25-Aug-2018.mat';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
load(file)

% load corresponding entropy file
entrFile = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/decoding/entropyDistr_N50M12_sp12_sig2_spW_25-Aug-2018.mat';
load(entrFile)
figure
errorbar(spAll',summH(:,1),summH(:,2),'o-','MarkerSize',12,'MarkerFaceColor',...
    Bu(10,:),'Color',Bu(10,:),'LineWidth',2,'CapSize',0)
lg = legend('N=50,M=12,sp=2,\sigma_c=2');
legend boxoff
set(lg,'FontSize',16)
xlabel('sparsity of W')
ylabel('differential entropy')

% meanClassAcur = mean(exp(-allCrossEntr),2);
% stdClassAcur = std(exp(-allCrossEntr),0,2);
nOdor = 50;
nRecp = 12;
noisSig = 0.1;
H = 500;    %hidden layer size
sp = 2;  
group = 2;  % groups of the data
%plot the figure
meanClassAcur = mean(summData.errorRate,2);
stdClassAcur = std(summData.errorRate,0,2);

figure
errorbar(summData.allSp',1-meanClassAcur,stdClassAcur,'o-','MarkerSize',10,'MarkerFaceColor',...
    Bu(9,:),'Color',Bu(9,:),'LineWidth',1.5,'CapSize',0)
lg = legend(['N=',num2str(nOdor),',M=',num2str(nRecp),',H=',num2str(H),',\Delta S=',num2str(noisSig)]);
set(lg,'FontSize',16)
ylim([0.7,1])
set(gca,'Xtick',0:0.2:1)
xlabel('$\rho_w$','interpreter','latex')
ylabel('accuracy')
figNamePref = ['class_Error_N',num2str(nOdor),'M',num2str(nRecp),...
    'sp',num2str(sp),'_nSig',num2str(noisSig),'_group',num2str(group),'_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% plot how classification accuracy changes with entropy
figure
scatter(summH(2:end,1),1-meanClassAcur,80,'o','LineWidth',2)
xlabel('differential entropy')
ylabel('classification accuracy')



% ================================================================
% Classification with different groups
% ================================================================
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding/class_diffGroup';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dFolder,filesep,'*.mat'));
filesRaw = {allFile.name}';

N = 100;
M = 12;
sp = 2;
sig = 2;
H = 500; %hidden layer size
nSig = 0.1;  %noise std
% screen the files, only part are used. Since the folder may contain
% heterogeneous data set
str_marker = 'classify_LDA';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

allGroup = 2:1:6;
allSp = 0.1:0.05:1;
meanError = zeros(length(allSp),length(allGroup));
stdError = zeros(length(allSp),length(allGroup));
minError = zeros(length(allGroup),2); % position and value
bestPerf = zeros(length(allGroup),2);  %store information of the optimal sparsity

for i0 = 1:length(files)
    s1 = '(?<= *group)[\d]+(?=_)';
    group = str2num(char(regexp(files{i0},s1,'match')));
    inx = find(allGroup == group);
    
    load(char(fullfile(dFolder,filesep,files{i0})));
    
    meanError(:,inx) = mean(summData.errorRate,2);
    stdError(:,inx) = std(summData.errorRate,0,2);
    
    [minX,minY] = sort(meanError(:,inx));
    minError(i0,2) = minX(1);  %minimum error
    minError(i0,1) = minY(1);  %minimum position
end

% plot the figure
figure
hold on
for i0 = 1:length(allGroup)
    plot(allSp',meanError(:,i0),'Color',Bu(2*i0 - 1,:))
    [ix, iy] = sort(meanError(:,i0));
    bestPerf(i0,1) = summData.allSp(iy(1));
    bestPerf(i0,2) = ix(1);
end
hold off
lg = legend('group = 2','group = 3','group = 4','group = 5','group = 6');
set(lg,'FontSize',16)
xlabel('$\rho_w$','interpreter','latex')
ylabel('Error')
figNamePref = ['class_Error_N',num2str(N),'M',num2str(M),...
    'sp',num2str(sp),'_nSig',num2str(sig),'_groupDpd_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])



% plot how the best performance change with group numbers
figure
plot(allGroup',minError(:,2),'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:))
xlabel('Number of groups')
ylabel('minimum error')
figNamePref = ['class_MinmumError_N',num2str(N),'M',num2str(M),...
    'sp',num2str(sp),'_nSig',num2str(sig),'_groupDpd_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% plot the optimal sparsity to achieve best performance
figure
plot(allGroup',bestPerf(:,1),'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:))
ylim([0,1])
xlabel('groups')
ylabel('$\rho_w^*$','interpreter','latex')
figNamePref = ['class_MinErrorPosi_N',num2str(N),'M',num2str(M),...
    'sp',num2str(sp),'_nSig',num2str(sig),'_groupDpd_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% ================================================================
% Classification performance with different receptors
% ================================================================
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/decoding';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
dName = 'classify_LDA_N100diffM_HS500_sp3_group2_npattern200_10-Aug-2018.mat';
load(fullfile(dFolder,dName))

allM = [5,8,10,12,15,20,25,30];
N = 100;
sp = 3;
sig = 2;

meanError = squeeze(mean(summData.errorRate,2));
stdError = squeeze(std(summData.errorRate,0,2));
bestPerf = zeros(size(meanError,2),2);  % best performance

figure
hold on
for i0 = 1:size(meanError,2)
    plot(summData.allSp',meanError(:,i0),'Color',seqColor(4 + 2*i0,:))
    [ix, iy] = sort(meanError(:,i0));
    bestPerf(i0,1) = summData.allSp(iy(1));
    bestPerf(i0,2) = ix(1);
end
lg = legend('M=5','M=8','M=10','M=12','M=15','M=20','M=25','M=30','Location','eastoutside');
set(lg,'FontSize',18)
legend boxoff

hold off
box on
xlabel('$\rho_w$','interpreter','latex')
ylabel('classification error')
figNamePref = ['class_Error_N',num2str(N),...
    'sp',num2str(sp),'_nSig',num2str(sig),'_diffM_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% best peformance, minimum error
figure
plot(allM',bestPerf(:,2),'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:))
box on
xlabel('$M$','interpreter','latex')
ylabel('minimum error')
figNamePref = ['class_MinmumError_N',num2str(N),...
    'sp',num2str(sp),'_nSig',num2str(sig),'_MDpd_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% sparsity that achieve best performance
figure
plot(allM',bestPerf(:,1),'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:))
box on
ylim([0,1])
xlabel('$M$','interpreter','latex')
ylabel('$\rho_w^*$','interpreter','latex')
figNamePref = ['class_MinErrorPosi_N',num2str(N),...
    'sp',num2str(sp),'_nSig',num2str(sig),'_MDpd_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])



% ================================================================
% Classification performance with different number of patterns
% ================================================================
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/decoding';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
dName = 'classify_LDA_N100M20diffPattern_HS500_sp3_group2_12-Aug-2018.mat';
load(fullfile(dFolder,dName))

allPattern = [10,20,50,100,200,500,1000];
N = 100;
M = 20;
sp = 3;
sig = 2;

meanError = squeeze(mean(summData.errorRate,2));
stdError = squeeze(std(summData.errorRate,0,2));
bestPerf = zeros(size(meanError,2),2);  % best performance

% the first data point seems to be wrong
slInx = 2:1:7;
figure
hold on
for i0 = 1:length(slInx)
    plot(summData.allSp',meanError(:,slInx(i0)),'Color',seqColor(4 + 2*slInx(i0),:))
    [ix, iy] = sort(meanError(:,slInx(i0)));
    bestPerf(slInx(i0),1) = summData.allSp(iy(1));
    bestPerf(slInx(i0),2) = ix(1);
end
lg = legend('cluster=20','cluster=50','cluster=100','cluster=200',...
    'cluster=500','cluster=1000','Location','eastoutside');
set(lg,'FontSize',18)
legend boxoff

hold off
box on
xlabel('$\rho_w$','interpreter','latex')
ylabel('classification error')
figNamePref = ['class_Error_N',num2str(N),...
    'sp',num2str(sp),'_nSig',num2str(sig),'_diffPattern_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% best peformance, minimum error
figure
plot(allPattern(slInx)',bestPerf(slInx,2),'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:))
box on
xlabel('clusters')
ylabel('minimum error')
figNamePref = ['class_MinError_N',num2str(N),...
    'sp',num2str(sp),'_nSig',num2str(sig),'_DiffPattern_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% sparsity that achieve best performance
figure
plot(allPattern(slInx)',bestPerf(slInx,1),'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:))
box on
ylim([0,1])
xlabel('clusters')
ylabel('$\rho_w^*$','interpreter','latex')
figNamePref = ['class_MinErrorPosi_N',num2str(N),...
    'sp',num2str(sp),'_nSig',num2str(sig),'_diffPattern_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% Classification with inhibition, comparison with excitation only situation
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
fName = 'LDA_inhi_r0__N50M10H200_sp2_G2_ns0.05_25-Dec-2018.mat';
refFile = 'LDA_N50M10sp2_ns0.05_nType2_P200_H200_dS0.1_25-Dec-2018.mat'; % reference
load(fullfile(dFolder,fName))

% load the excitation-only data
load(fullfile(dFolder,refFile))
aveError = mean(allTestError,2);
stdError = std(allTestError,0,2);
aveMin = min(aveError);



nOdor = 50;
nRecp = 10;
noisSig = 0.05;
dS = 0.1;
H = 200;
NP = 500;

figure
hold on
errorbar(summData.allSp',mean(summData.errorRate,2),std(summData.errorRate,0,2),'o-','MarkerSize',10,'MarkerFaceColor',...
    Bu(9,:),'Color',Bu(9,:),'LineWidth',1.5,'CapSize',0)
% set(lg,'FontSize',16)
plot([0;1],[aveMin;aveMin],'k--')
ylim([0.05,0.2])
box on
set(gca,'Xtick',0:0.2:1)
xlabel('$r_0$','interpreter','latex')
ylabel('classification error')
% figNamePref = ['class_Error_N',num2str(nOdor),'M',num2str(nRecp),...
%     'sp',num2str(sp),'_nSig',num2str(noisSig),'_group',num2str(group),'_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

figure
errorbar(allSp',aveError,stdError,'o-','MarkerSize',10,'MarkerFaceColor',...
    Bu(9,:),'Color',Bu(9,:),'LineWidth',1.5,'CapSize',0)
xlabel('$\rho_w$','interpreter','latex')
ylabel('classification error')

%% classification, SVM method
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
dName = 'classify_svm_N100M10H500_sp3_G2_ns0.05_24-Sep-2018.mat';
load(fullfile(dFolder,dName))

N = 100;
M = 10;
sp = 3;
sig = 2;
noiseSig = 0.05;

allSp = 0.05:0.05:0.95;


figure
Bu = brewermap(11,'Blues');    % blues
errorbar(allSp',mean(summData.errorRate,2),std(summData.errorRate,0,2),...
        'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:),'LineWidth',2,'CapSize',0)
lg = legend('N = 100,M=10,\sigma_n = 0.05, P = 100');
legend boxoff
set(lg, 'FontSize',16)

xlabel('$\rho_w^*$','interpreter','latex')
ylabel('MSE')
figNamePref = ['class_SVM_MinErrorPosi_N',num2str(N),'M',num2str(M),...
    'sp',num2str(sp),'_nSig',num2str(noiseSig),'_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% plot the alpha - m when N = 2, M goes infinity
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
file = 'alpha_m.fig';

openfig(fullfile(dFolder,file));
fh = gcf;
% ah = fh.Children;
axesObjs = get(fh, 'Children');        %axes handles
dataObjs = get(axesObjs, 'Children');  %handles to low-level graphics objects in axes
data = dataObjs(2); 
lineData = cell(3,1);
ScatterData = cell(3,1);
for i0 = 1:3
    temp = data{:};
    lineData{i0}(:,1) = temp(i0).XData;
    lineData{i0}(:,2) = temp(i0).YData;
    ScatterData{i0}(:,1) = temp(3 + i0).XData;
    ScatterData{i0}(:,2) = temp(3 + i0).YData;
end

% replot the figure
figure
hold on
for i0 = 1:3
    plot(ScatterData{i0}(:,1),ScatterData{i0}(:,2),'o-','Color',Bu(2 + 3*i0,:),...
        'MarkerFaceColor',Bu(2 + 3*i0,:),'MarkerSize',8,'LineWidth',1)
end
box on

lg = legend('\sigma_c = 1','\sigma_c = 2','\sigma_c = 3');
set(lg,'FontSize',18)
legend boxoff
for i0 = 1:3
    plot(lineData{i0}(:,1),lineData{i0}(:,2),'k--','lineWidth',1.5)
end
hold off

set(gca,'XScale','log')
xlabel('$\ln(N)$','Interpreter','latex')
ylabel('$\alpha$','Interpreter','latex')

figNamePref = ['N2M_alpha_m_sigc_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% Gcmi with correlated concentration input, N dependent
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% load the summary data set
dName = 'gcmiCorr_diffN_M10sp3sig2_13-Aug-2018.mat';
load(fullfile(dFolder,dName))

M = 10;
sp = 3;
sig = 2;

% ================================
% fmin changes with N
% ================================
figure
errorbar(allN,allMeanFmin,allStdFmin,'o-','MarkerSize',12,'Color',Bu(9,:), ...
    'MarkerFaceColor',Bu(9,:),'LineWidth',2,'CapSize',0)
lg = legend(['M=',num2str(M),',sp=',num2str(sp),',\sigma_c =',num2str(sig)],'Location','southeast');
set(lg,'FontSize',20)
legend boxoff
xlabel('$N$','interpreter','latex')
ylabel('$f_{min}$','interpreter','latex')
prefix = ['GcmiCorr_fmin_Ndp_M',num2str(M),'sig',num2str(sig),'_',date,'.mat'];
print(gcf,'-depsc',fullfile(sFolder,[prefix,'.eps']))
saveas(gcf,fullfile(sFolder,[prefix,'.fig']))

% ================================
% sp W changes with N
% ================================
figure
errorbar(allN,allMeanSpW,allStdSpW,'o-','MarkerSize',12,'Color',Bu(9,:), ...
    'MarkerFaceColor',Bu(9,:),'LineWidth',2,'CapSize',0)
lg = legend(['M=',num2str(M),',sp=',num2str(sp),',\sigma_c =',num2str(sig)],'Location','southeast');
set(lg,'FontSize',20)
legend boxoff
xlabel('$N$','interpreter','latex')
ylabel('$\rho_{w}$','interpreter','latex')
prefix = ['GcmiCorr_spW_Ndp_M',num2str(M),'sig',num2str(sig),'_',date,'.mat'];
print(gcf,'-depsc',fullfile(sFolder,[prefix,'.eps']))
saveas(gcf,fullfile(sFolder,[prefix,'.fig']))

% ================================
% mu W changes with N
% ================================
figure
errorbar(allN,allMeanAveW,allStdAveW,'o-','MarkerSize',12,'Color',Bu(9,:), ...
    'MarkerFaceColor',Bu(9,:),'LineWidth',2,'CapSize',0)
lg = legend(['M=',num2str(M),',sp=',num2str(sp),',\sigma_c =',num2str(sig)],'Location','southeast');
set(lg,'FontSize',20)
ylim([-2,-1])
legend boxoff
xlabel('$N$','interpreter','latex')
ylabel('$\mu_{w}$','interpreter','latex')
prefix = ['GcmiCorr_meanW_Ndp_M',num2str(M),'sig',num2str(sig),'_',date,'.mat'];
print(gcf,'-depsc',fullfile(sFolder,[prefix,'.eps']))
saveas(gcf,fullfile(sFolder,[prefix,'.fig']))

% ================================
% std W changes with N
% ================================
figure
errorbar(allN,allMeanSigW,allStdSigW,'o-','MarkerSize',12,'Color',Bu(9,:), ...
    'MarkerFaceColor',Bu(9,:),'LineWidth',2,'CapSize',0)
lg = legend(['M=',num2str(M),',sp=',num2str(sp),',\sigma_c =',num2str(sig)],'Location','southeast');
set(lg,'FontSize',20)
legend boxoff
ylim([1,2])
xlabel('$N$','interpreter','latex')
ylabel('$\sigma_{w}$','interpreter','latex')
prefix = ['GcmiCorr_meanW_Ndp_M',num2str(M),'sig',num2str(sig),'_',date,'.mat'];
print(gcf,'-depsc',fullfile(sFolder,[prefix,'.eps']))
saveas(gcf,fullfile(sFolder,[prefix,'.fig']))

%% Debug with GCMI, different initial w0
% this part compare the sigmac-depenence of optimal W when different
% initial conditions were used
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% load the summary data set
dName = 'N100M30sp2Summary_h1_2018-08-10.mat';
load(fullfile(dFolder,dName))

% extra parameters
N = 100;
M = 30;
sp = 2;
sig = 2;

figure
subplot(2,2,1)
errorbar(dataSumm.allSig',dataSumm.spMean,dataSumm.spStd,'o-','MarkerSize',12,'Color',Bu(9,:), ...
    'MarkerFaceColor',Bu(9,:),'LineWidth',2,'CapSize',0)
xlabel('$\sigma_c$','interpreter','latex')
ylabel('$\rho_w$','interpreter','latex')


subplot(2,2,2)
errorbar(dataSumm.allSig',dataSumm.meanAveW,dataSumm.stdAveW,'o-','MarkerSize',12,'Color',Bu(9,:), ...
    'MarkerFaceColor',Bu(9,:),'LineWidth',2,'CapSize',0)
xlabel('$\sigma_c$','interpreter','latex')
ylabel('$\mu_w$','interpreter','latex')


subplot(2,2,3)
errorbar(dataSumm.allSig',dataSumm.meanStdW,dataSumm.stdStdW,'o-','MarkerSize',12,'Color',Bu(9,:), ...
    'MarkerFaceColor',Bu(9,:),'LineWidth',2,'CapSize',0)
xlabel('$\sigma_c$','interpreter','latex')
ylabel('$\sigma_w$','interpreter','latex')


subplot(2,2,4)
errorbar(dataSumm.allSig',dataSumm.allFmean,dataSumm.allFstd,'o-','MarkerSize',12,'Color',Bu(9,:), ...
    'MarkerFaceColor',Bu(9,:),'LineWidth',2,'CapSize',0)
xlabel('$\sigma_c$','interpreter','latex')
ylabel('$f_{min}$','interpreter','latex')


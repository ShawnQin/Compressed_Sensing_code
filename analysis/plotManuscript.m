  % this program plots the figures in manuscripts
% we load the summarized data and combine them to form a concret
% representation
% this is modified from finalPlots.m

% last revised on 07/19/2018

close all
clear

%% set some color and default graphical settings
% this section should be ran before the following sections
dFolder = '/Users/shan/Dropbox/olfactionProject/data/gcmi';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% graphics settings
defaultGraphicsSetttings

%define some colors using brewermap
RdBu = brewermap(11,'RdBu');   % red and blue
Bu = brewermap(11,'Blues');    % blues
Gr = brewermap(11,'Greys');    % greys

myBu = [3, 110, 184]/256;
%myOr = [224, 135, 51]/256;
myOr = [255, 185, 85]/256;

myRd = [202, 39, 44]/256;
lBu = [96,166,223]/255; %light blue
dpBu = [63,114,183]/255; % deep blue
dkBu = [50,78,147]/255;   %dark blue
Or = [220,150,71]/255;  % orange
brickRd = [201,69,89]/255;  %brick red
green = [107,169,81]/255;  %green
purple = [113,49,119]/255;  % purple

%% Figure 1C nonlinear ORN
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

C = 10.^(-2:0.01:2);
Km = 1;
n = 1;
f = @(Km,C) (C/Km).^n./(1+(C/Km).^n);

figure
plot(C,f(Km,C),'k','LineWidth',4)
set(gca,'xscale','log','XTick',10.^(-2:1:2))
xlabel('$\log_{10}(c)$','Interpreter','latex')
ylabel('normalized repsonse')

prefix = ['nonlinearORN_h1',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% Figure 3 how sparsity of W change with sigma_c, n and N
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% ===============================================================
% define the layout of the graphics, we have 3x3 subplots
% ===============================================================
errorBarSize = 1.5;
errorMarkerSize  = 10;
LineWidth = 1.5;
labelFontSize = 24;
axisFontSize = 20;
ticketWidth = 1.5;

colorInx = 1;   % the color index in a 11 order of colors

% figureSize = [0 0 9 7];
figureSize = [0 0 11 7];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(3,3,[.03 .08],[.1 .04],[.1 .05]);


% =========================================================
% summarized data for different sigma, we only need N = 100
% =========================================================

% fName = 'gcmi_diffN_summary_06222018.mat';
% fName = 'gcmi_Ndp_M20Sp2sig2_01-Sep-2018.mat';
% load(fullfile(dFolder,fName));

%{
Ninx = 5;  %the index corresponding to N = 100

% sp of W vs sigmac
axes(ha(1))
errorbar(allSig,meanSp(:,Ninx),stdSp(:,Ninx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('N = 100,n = 2','Location','southeast');
legend boxoff
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

% mean W vs sigmac
axes(ha(4))
errorbar(allSig,meanPeak(:,Ninx),stdPeak(:,Ninx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

%std w vs sigmac
axes(ha(7))
errorbar(allSig,meanSig(:,Ninx),stdSig(:,Ninx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_w$','Interpreter','latex','FontSize',labelFontSize)
% ylabel('mean W','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')
%}

% =========================================================
% summarized data for different sigma, we only need N = 100
% =========================================================

% fName = 'N100M20sp2_Summary_diffSig_2018-08-25.mat';
% fName = 'N50M13sp2_Summary_diffSig_2018-10-05.mat';
fName = 'gcmi_sigdp_N50M13Sp3_19-Oct-2018.mat';
load(fullfile(dFolder,fName));

inx = 4:1:17;  %the index corresponding to N = 100

% sp of W vs sigmac
axes(ha(1))
errorbar(dataSumm.allSig(inx)',dataSumm.meanSpW(inx),dataSumm.stdSpW(inx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('N = 50,n = 3','Location','southeast');
legend boxoff
% ylim([0.6,0.8])
% xlim([1,4.1])
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

% mean W vs sigmac
axes(ha(4))
errorbar(dataSumm.allSig(inx)',dataSumm.meanW(inx),dataSumm.stdW(inx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% xlim([1,4.1])
% ylim([-3.5,0.5])
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

%std w vs sigmac
axes(ha(7))
errorbar(dataSumm.allSig(inx)',dataSumm.meanSigW(inx),dataSumm.stdSigW(inx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
xlim([2,6])
% ylim([1,2.2])
xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_w$','Interpreter','latex','FontSize',labelFontSize)
% ylabel('mean W','FontSize',labelFontSize)
set(gca,'XTick',2:1:6,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')


% ======================================================
% summarized data for different input sparsity (n) 
% ======================================================
% nFile = 'N100M20Summary_h1_28-Jun-2018.mat';
% nFile = 'gcmi_Spdpd_N100M20Summary_01-Sep-2018.mat';
nFile = 'gcmi_spWdp_N50M13Sp9_19-Oct-2018.mat';
load(fullfile(dFolder,nFile))

N = 50;
M = 13;
sig = 2;

inx = 1:1:5;

axes(ha(2))
errorbar(dataSumm.sp(inx)',dataSumm.meanSpW(inx),dataSumm.stdSpW(inx),'o-','MarkerSize',...
      errorMarkerSize,'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'CapSize',0,'LineWidth',LineWidth)
lg = legend(['N = 50,\sigma_c = ',num2str(sig)],'Location','southwest');
legend boxoff
% ylim([0.4,0.9])
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

axes(ha(5))
errorbar(dataSumm.sp(inx)',dataSumm.meanW(inx),dataSumm.stdW(inx),'o-','MarkerSize',...
      errorMarkerSize,'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),...
      'CapSize',0,'LineWidth',LineWidth)
% ylim([-3.5,0.5])
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

axes(ha(8))
errorbar(dataSumm.sp(inx)',dataSumm.meanSigW(inx),dataSumm.stdSigW(inx),'o-','MarkerSize',...
      errorMarkerSize,'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'CapSize',0,'LineWidth',LineWidth)
ylim([1,1.5])
xlabel('$n$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_w$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',2:1:6)



% ======================================================
% summarized data for different number of odorants (N) 
% ======================================================
% nFile = 'gcmi_Ndp_M20Sp2sig2_01-Sep-2018.mat';
nFile = 'gcmi_Ndp_M13Sp3sig2_19-Oct-2018.mat';

%nFile = 'gcmi_Ndp_M20Sp2sig2_12-Sep-2018.mat';

load(fullfile(dFolder,nFile))

M = 13;
sp = 3;
sig = 2;
allN = [40:10:100,120];

odorInx = 1:1:8;   %odor selected
% odorInx = 1:1:8; %odor selected

% sp of W vs sigmac
axes(ha(3))
errorbar(dataSumm.allN(odorInx)',dataSumm.meanSpW(odorInx),dataSumm.stdSpW,'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('n=2,\sigma_c = 2','Location','southwest');
legend boxoff
% ylim([0.4,0.9])
% xlim([50,150])
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

% mean W vs sigmac
axes(ha(6))
errorbar(dataSumm.allN(odorInx)',dataSumm.meanW(odorInx),dataSumm.stdW(odorInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
ylim([-1.5,-1])
xlim([40 120])
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

%std w vs sigmac
axes(ha(9))
errorbar(dataSumm.allN(odorInx)',dataSumm.meanSigW(odorInx),dataSumm.stdSigW(odorInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% ylim([1.2,1.7])
% ylim([1,2.2])
xlim([40 120])
xlabel('$N$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_w$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

%save the figure
prefix = ['gcmi_summ_figure3_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% Figure 3 how sparsity of W change with sigma_c, n and N, reorganized
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% ===============================================================
% define the layout of the graphics, we have 3x3 subplots
% ===============================================================
errorBarSize = 1.5;
errorMarkerSize  = 10;
LineWidth = 1.5;
labelFontSize = 24;
axisFontSize = 20;
ticketWidth = 1.5;

colorInx = 1;   % the color index in a 11 order of colors

% figureSize = [0 0 9 7];
figureSize = [0 0 5.1 3.4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(2,1,[.0 .2],[0.15,0.03],[0.15,0.03]);


% =========================================================
% summarized data for different sigma, we only need N = 100
% =========================================================

fName = 'gcmi_sigdp_N50M13Sp3_19-Oct-2018.mat';
load(fullfile(dFolder,fName));

inx = 4:1:17;  %the index corresponding to N = 100

% sp of W vs sigmac
axes(ha(1))
errorbar(dataSumm.allSig(inx)',dataSumm.meanSpW(inx),dataSumm.stdSpW(inx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('N = 50,n = 3','Location','southeast');
legend boxoff
% ylim([0.6,0.8])
% xlim([1,4.1])
xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% mean W vs sigmac
axes(ha(1))
errorbar(dataSumm.allSig(inx)',dataSumm.meanW(inx),dataSumm.stdW(inx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% xlim([1,4.1])
% ylim([-3.5,0.5])
% xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

%std w vs sigmac
axes(ha(2))
errorbar(dataSumm.allSig(inx)',dataSumm.meanSigW(inx),dataSumm.stdSigW(inx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
xlim([2,6])
% ylim([1,2.2])
xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_w$','Interpreter','latex','FontSize',labelFontSize)
% ylabel('mean W','FontSize',labelFontSize)
set(gca,'XTick',2:1:6,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')


% ======================================================
% summarized data for different input sparsity (n) and N
% ======================================================
figureSize = [0 0 8 4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(1,2,[.12 .1],[.2 .01],[.1 .02]);

nFile = 'gcmi_spWdp_N50M13Sp9_19-Oct-2018.mat';
load(fullfile(dFolder,nFile))

N = 50;
M = 13;
sig = 2;

inx = 1:1:5;

axes(ha(1))
errorbar(dataSumm.sp(inx)',dataSumm.meanSpW(inx),dataSumm.stdSpW(inx),'o-','MarkerSize',...
      errorMarkerSize,'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'CapSize',0,'LineWidth',LineWidth)
lg = legend(['N = 50,\sigma_c = ',num2str(sig)],'Location','southwest');
legend boxoff
% ylim([0.4,0.9])
xlabel('$n$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')


% ======================================================
% summarized data for different number of odorants (N) 
% ======================================================
nFile = 'gcmi_Ndp_M13Sp3sig2_19-Oct-2018.mat';
load(fullfile(dFolder,nFile))

M = 13;
sp = 3;
sig = 2;
allN = [40:10:100,120];

odorInx = 1:1:8;   %odor selected

axes(ha(2))
errorbar(dataSumm.allN(odorInx)',dataSumm.meanSpW(odorInx),dataSumm.stdSpW,'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('n=2,\sigma_c = 2','Location','southwest');
legend boxoff
xlabel('$N$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

%save the figure
prefix = ['gcmi_summ_figure3_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% Figure 3 summary, different layout 1
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% ===============================================================
% define the layout of the graphics, we have 3x3 subplots
% ===============================================================
errorBarSize = 1.5;
errorMarkerSize = 10;
LineWidth = 1.5;
labelFontSize = 24;
axisFontSize = 20;
ticketWidth = 1.5;

colorInx = 1;   % the color index in a 11 order of colors

figureSize = [0 0 4.5 3.4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
fName = 'gcmi_sigdp_N50M13Sp3_19-Oct-2018.mat';
load(fullfile(dFolder,fName));

inx = 4:1:17;  %the index corresponding to N = 100

[ha, pos] = tight_subplot(2,1,[.0 .2],[0.15,0.03],[0.15,0.03]);
axes(ha(1))
errorbar(dataSumm.allSig(inx)',dataSumm.meanW(inx),dataSumm.stdW(inx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

%std w vs sigmac
axes(ha(2))
errorbar(dataSumm.allSig(inx)',dataSumm.meanSigW(inx),dataSumm.stdSigW(inx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
xlim([2,6])
% ylim([1,2.2])
xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_w$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'XTick',2:1:6,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')


% figure of how rho_w changes with sigma_c
figureSize = [0 0 4.5 3.4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
errorbar(dataSumm.allSig(inx)',dataSumm.meanSpW(inx),dataSumm.stdSpW(inx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('N = 50,n = 3','Location','southeast');
legend boxoff
% ylim([0.6,0.8])
% xlim([1,4.1])
xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% ============================================================
% how rho_w changes with input sparsity
% ============================================================
nFile = 'gcmi_spWdp_N50M13Sp9_19-Oct-2018.mat';
load(fullfile(dFolder,nFile))
figureSize = [0 0 4.5 3.4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
inx = 1:1:5;
errorbar(dataSumm.sp(inx)',dataSumm.meanSpW(inx),dataSumm.stdSpW(inx),'o-','MarkerSize',...
      errorMarkerSize,'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'CapSize',0,'LineWidth',LineWidth)
lg = legend(['N = 50,\sigma_c = ',num2str(sig)],'Location','southwest');
legend boxoff
% ylim([0.4,0.9])
xlabel('$n$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')


% ============================================================
% how rho_w changes with total number of odorants
% ============================================================
nFile = 'gcmi_Ndp_M13Sp3sig2_19-Oct-2018.mat';
load(fullfile(dFolder,nFile))
figureSize = [0 0 4.5 3.4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');


allN = [40:10:100,120];
odorInx = 1:1:8;   %odor selected
errorbar(dataSumm.allN(odorInx)',dataSumm.meanSpW(odorInx),dataSumm.stdSpW,'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('n=3,\sigma_c = 2','Location','southwest');
legend boxoff
xlabel('$N$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% ============================================================
% how rho_w changes with number of receptors
% ============================================================
% nFile = 'gcmi_Mdp_N100Sp3sig2_14-Dec-2018.mat';
nFile = 'gcmi_Mdp_N50Sp3sig2_15-Dec-2018.mat';
load(fullfile(dFolder,nFile))
figureSize = [0 0 4.5 3.4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

N = 50;
% allM = 5:5:40;
allM = [10,20:5:50];
orInx = 2:1:8;   %odor selected

figure
errorbar(dataSumm.allM(orInx)',dataSumm.meanSpW(orInx),dataSumm.stdSpW(orInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,'CapSize',0)
lg = legend('N = 100,n=3,\sigma_c = 2','Location','southwest');
legend boxoff
ylim([0.5,0.62])
xlabel('$M$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')
prefix = ['gcmi_fig3_diffM_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% how mu_w and sigma_w change with input
figure
orInx = 2:1:8;   %odor selected
errorbar(dataSumm.allM(orInx)',dataSumm.meanW(orInx),dataSumm.stdW(orInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,'CapSize',0)
ylim([-1.4,-1])
xlabel('$M$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)

figure
orInx = 2:1:8;   %odor selected
errorbar(dataSumm.allM(orInx)',dataSumm.meanSigW(orInx),dataSumm.stdSigW(orInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,'CapSize',0)
ylim([0.9,1.3])
xlabel('$M$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)

%% Figure 3 summary, different layout  2
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% ===============================================================
% define the layout of the graphics, we have 3x3 subplots
% ===============================================================
errorBarSize = 1.5;
errorMarkerSize  = 10;
LineWidth = 1.5;
labelFontSize = 24;
axisFontSize = 20;
ticketWidth = 1.5;

% basic parameters
N = 50;
M = 13;
sp = 3;
sig = 2;

colorInx = 1;   % the color index in a 11 order of colors

% sigma-dependent
fName = 'gcmi_sigdp_N50M13Sp3_19-Oct-2018.mat';
load(fullfile(dFolder,fName));
summSigdp = dataSumm;

% N-dependent
% NFile = 'gcmi_Ndp_M13Sp3sig2_19-Oct-2018.mat';
NFile = 'gcmi_Ndp_M13Sp3sig2_24-Oct-2018.mat';
load(fullfile(dFolder,NFile))
summNdp = dataSumm;

% sp-dependent
nFile = 'gcmi_spWdp_N50M13Sp9_19-Oct-2018.mat';
load(fullfile(dFolder,nFile))
summSpdp = dataSumm;

% M-dependent
Mfile = 'gcmi_Mdp_N50Sp3sig2_15-Dec-2018.mat';
load(fullfile(dFolder,Mfile))
summMdp = dataSumm;

% ========================================================
% first row plot how rho_w change with sigma_c, n, N, M
% ========================================================
figureSize = [0 0 16 3];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
[ha, pos] = tight_subplot(1,4,[.01 .08],[.23 .04],[.08 .01]);


inx1 = 4:1:17;  %the index corresponding to N = 100

axes(ha(1))
errorbar(summSigdp.allSig(inx1)',summSigdp.meanSpW(inx1),summSigdp.stdSpW(inx1),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend(['n = ',num2str(sp)],'Location','northwest');
legend boxoff
xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_w$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

axes(ha(2))
inx2 = 1:1:5;
errorbar(summSpdp.sp(inx2)', summSpdp.meanSpW(inx2),summSpdp.stdSpW(inx2),'o-','MarkerSize',...
      errorMarkerSize,'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'CapSize',0,'LineWidth',LineWidth)
lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
legend boxoff
% ylim([0.4,0.9])
xlabel('$n$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% allN = [40,50,70,90,100,120];
odorInx = [1,2,4:1:8];   %odor selected
axes(ha(3))
errorbar(summNdp.allN(odorInx)',summNdp.meanSpW(odorInx),summNdp.stdSpW(odorInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('n=3,\sigma_c = 2','Location','northeast');
legend boxoff
ylim([0.57,0.72])
xlabel('$N$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

MInx = 2:1:8;
axes(ha(4))
errorbar(summMdp.allM(MInx)',summMdp.meanSpW(MInx),summMdp.stdSpW(MInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',purple,'Color',purple,'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('N = 50,n=3,\sigma_c = 2','Location','northeast');
legend boxoff
ylim([0.5,0.62])
xlabel('$M$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

prefix = ['gcmi_summ_figure3_rhow_all',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])



% ========================================================
% second row plot how sigma_w change with sigma_c, n, N, M
% ========================================================
figureSize = [0 0 16 3];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
[ha, pos] = tight_subplot(1,4,[.01 .08],[.23 .04],[.08 .01]);


inx1 = 4:1:17;  %the index corresponding to N = 100

axes(ha(1))
errorbar(summSigdp.allSig(inx1)',summSigdp.meanSigW(inx1),summSigdp.stdSigW(inx1),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend(['n = ',num2str(sp)],'Location','northwest');
legend boxoff
xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_w$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

axes(ha(2))
inx2 = 1:1:5;
errorbar(summSpdp.sp(inx2)', summSpdp.meanSigW(inx2),summSpdp.stdSigW(inx2),'o-','MarkerSize',...
      errorMarkerSize,'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'CapSize',0,'LineWidth',LineWidth)
lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
legend boxoff
ylim([1,1.5])
xlabel('$n$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% allN = [40,50,70,90,100,120];
odorInx = [1,2,4:1:8];   %odor selected
axes(ha(3))
errorbar(summNdp.allN(odorInx)',summNdp.meanSigW(odorInx),summNdp.stdSigW(odorInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('n=3,\sigma_c = 2','Location','northeast');
legend boxoff
ylim([0.9,1.55])
xlabel('$N$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

MInx = 2:1:8;
axes(ha(4))
errorbar(summMdp.allM(MInx)',summMdp.meanSigW(MInx),summMdp.stdSigW(MInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',purple,'Color',purple,'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('N = 50,n=3,\sigma_c = 2','Location','northeast');
legend boxoff
ylim([0.9,1.4])
xlabel('$M$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

prefix = ['gcmi_summ_figure3_sigW_all',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% ================================================================
% mu_w and rho_w, different sigma_c
% ================================================================
figureSize = [0 0 3.8 3];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(2,1,[.0 .2],[0.2,0.03],[0.2,0.07]);

axes(ha(1))
errorbar(summSigdp.allSig(inx1)',summSigdp.meanW(inx1),summSigdp.stdW(inx1),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
legend boxoff
% xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

axes(ha(2))
errorbar(summSigdp.allSig(inx1)',summSigdp.meanSigW(inx1),summSigdp.stdSigW(inx1),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
% legend boxoff
xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')


% ================================================================
% mu_w and rho_w, different n
% ================================================================
figureSize = [0 0 3.8 3];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(2,1,[.0 .2],[0.2,0.03],[0.2,0.07]);

axes(ha(1))
errorbar(summSpdp.sp(inx2)',summSpdp.meanW(inx2),summSpdp.stdW(inx2),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
% legend boxoff
% xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

axes(ha(2))
errorbar(summSpdp.sp(inx2)',summSpdp.meanSigW(inx2),summSpdp.stdSigW(inx2),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
% legend boxoff
ylim([1,1.5])
xlabel('$n$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')


% ================================================================
% mu_w and rho_w, different N
% ================================================================
figureSize = [0 0 3.8 3];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(2,1,[.0 .2],[0.2,0.03],[0.2,0.07]);

axes(ha(1))
errorbar(summNdp.allN(odorInx)',summNdp.meanW(odorInx),summNdp.stdW(odorInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
% legend boxoff
ylim([-1.6,-1])
% xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

axes(ha(2))
errorbar(summNdp.allN(odorInx)',summNdp.meanSigW(odorInx),summNdp.stdSigW(odorInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
% legend boxoff
% ylim([1,1.5])
xlabel('$N$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')


% ================================================================
% mu_w and rho_w, different M
% ================================================================
figureSize = [0 0 3.8 3];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(2,1,[.0 .2],[0.2,0.03],[0.2,0.07]);

MInx = 2:1:8;
axes(ha(1))
errorbar(summMdp.allM(MInx)',summMdp.meanW(MInx),summMdp.stdW(MInx),'o-',...
    'MarkerSize',errorMarkerSize,'MarkerFaceColor',purple,'Color',purple,...
    'LineWidth',LineWidth,'CapSize',0)

ylim([-1.4,-1])
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

axes(ha(2))
errorbar(summMdp.allM(MInx)',summMdp.meanSigW(MInx),summMdp.stdSigW(MInx),...
    'o-','MarkerSize',errorMarkerSize,'MarkerFaceColor',purple,...
    'Color',purple,'LineWidth',LineWidth,'CapSize',0)
ylim([0.9,1.3])
xlabel('$M$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

%% Figure 4 Mean Field Theory with N = 2, direct integration
dFolder = '/Users/shan/Dropbox/olfactionProject/data/tempfigure';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
dName = 'data_fit_MFT.mat'; %more direct and accuarate, see SI
% dName = 'MFT_fit.mat';  % Yuhai's ansatz
load(fullfile(dFolder, dName))

inx = [1 3 6];
allsig = 1.1:0.2:2.9; % all the sigma

[Y1,Y2] = sort(all,2,'descend');
largest = zeros(size(all,1),2);
for i0 = 1:size(largest,1)
    largest(i0,1) = q_all(Y2(i0,1));
end
largest(:,2) = Y1(:,1);
% ==================================================
% plot three I2(m) with different sigma_;c
% ==================================================
figure
hold on
for i0 = 1:length(inx)
    plot(q_all',all(inx(i0),:)','Color',Bu(2+i0*3,:),'LineWidth',3)
    ah = gca;
    yl = ah.YLim;
    plot([largest(inx(i0),1);largest(inx(i0),1)],yl','--','LineWidth',2,'Color',Bu(2+i0*3,:))
end
ah.YLim = [0 5];
lg = legend('\sigma_c = 1.1','\sigma_c = 1.5','\sigma_c = 2.1');
set(lg,'FontSize',16)
legend boxoff
box on

xlabel('$\rho_w$','interpreter','latex')
ylabel('$I_2$','interpreter','latex')
hold off
prefix = ['fig4_MFT_N2_entr_rhoW',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ==================================================
% plot peak position of I2 with differnt sigma_c
% ==================================================
figure
plot(allsig',largest(:,1),'o','MarkerSize',12,'MarkerFaceColor',myOr,...
    'MarkerEdgeColor','k','LineWidth',1)
ylim([0.6 0.75])
prefix = ['fig4_MFT_N2_rhow_sigc',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% Figure 4 Mean Field Theory with N = 2, Yuhai's ansatz
dFolder = '/Users/shan/Dropbox/olfactionProject/data/tempfigure';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
dName = 'MFT_fit.mat';  % Yuhai's ansatz
load(fullfile(dFolder, dName))

inx = [1 3 6];
allsig = 1.1:0.2:2.9; % all the sigma
allM = [60,100,150,200];
ixM = 3;   %select M = 150 to plot
allRho = 0.01:0.01:1;

figure
hold on
for i0 = 1:length(inx)
    plot(allRho,fitted_fcn{ixM,inx(i0)},'Color',Bu(2+i0*3,:),'LineWidth',3)
    ah = gca;
    yl = ah.YLim;
    plot([fitted_q(ixM,inx(i0));fitted_q(ixM,inx(i0))],yl','--','LineWidth',2,'Color',Bu(2+i0*3,:))
end
% ah.YLim = [2 4.5];
% xlim([0.1,1])
lg = legend('\sigma_c = 1.1','\sigma_c = 1.5','\sigma_c = 2.1');
set(lg,'FontSize',16)
legend boxoff
box on
xlabel('$\rho_w$','interpreter','latex')
ylabel('$I_2$','interpreter','latex')
hold off
prefix = ['fig4_MFT_N2_entr_rhoW',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ==================================================
% plot peak position of I2 with differnt sigma_c
% ==================================================
figure
plot(allsig,fitted_q(ixM,:),'o','MarkerSize',12,'MarkerFaceColor',myOr,...
    'MarkerEdgeColor','k','LineWidth',1)
ylim([0.6 0.95])
prefix = ['fig4_MFT_N2_rhow_sigc',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% Supplemental Figure Mean Field Theory
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% ===============================================================
% define the layout of the graphics, we have 3x3 subplots
% ===============================================================
errorMarkerSize  = 10;
LineWidth = 2;
labelFontSize = 28;
axisFontSize = 24;
ticketWidth = 1.5;

colorInx = 1;   % the color index in a 11 order of colors

figureSize = [0 0 13 4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(1,3,[.07 .07],[.25 .05],[.1 .03]);


fName = 'MFTsigmac_summ.mat';
load(fullfile(dFolder,fName));

% sp of W vs sigmac
axes(ha(1))
plot(sigall,sparsity,'Color',Bu(10,:))
% lg = legend('N = 100,n = 2','Location','southeast');
% legend boxoff
% xlim([1,5])
xlabel('$\sigma_{c}$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'XTick',1:1:5,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% mean W vs sigmac
axes(ha(2))
plot(sigall,meanvalue,'Color',Bu(10,:))
% xlim([1,5])
xlabel('$\sigma_{c}$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'XTick',1:1:5,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

%std w vs sigmac
axes(ha(3))
plot(sigall,stdvalue,'Color',Bu(10,:))
% xlim([1,5])
xlabel('$\sigma_{c}$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'XTick',1:1:5,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

%save the figure
prefix = ['MFT_figure5_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% Figure 5 Reconstruction performance
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% load the data
dName = 'recons_tensor_N100M20sp2sig2ns0.05_loss_14-Oct-2018New.mat';
load(fullfile(dFolder,dName))

% load plastic reconstruction, first layer is "learned"
temp = load('/Users/shan/Documents/machine learning/olfactionEnco/loss_plastic_100M20noise0.05_sp0.55_L2.mat');
plsLoss = [mean(temp.testLost),std(temp.testLost)];
% add extral simulation data, 20 repeats, lack the sp = 1
% temp = load(fullfile(dFolder,'recons_tensor_N100M20sp2sig2ns0.05_loss_03-Oct-2018.mat'));

% cmbData = nan(21,50);
% cmbData(2:end,:) = [allTestLoss(2:end,:), temp.allTestLoss];
% cmbData(1,1:40) = allTestLoss(1,:);
% 
entrData = 'entropyDistr_N100M20_sp2_sig2_ns0.05_03-Oct-2018.mat';
load(fullfile(dFolder,entrData))


% parameters
N = 100;
M = 20;
spar = 2;
sig = 2;
noise = 0.05;

% allSp = 1:-0.05:0;
inx  = 1:1:18;

figure
hold on
% errorbar(allSp(inx)',nanmean(log10(cmbData(inx,:)),2),nanstd(log10(cmbData(inx,:)),0,2),...
%         'o-','MarkerSize',10,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:),'LineWidth',1.5,'CapSize',0)
errorbar(1-allSp(inx)',nanmean(log10(allTestLoss(inx,:)),2),nanstd(log10(allTestLoss(inx,:)),0,2),...
        'o-','MarkerSize',10,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:),'LineWidth',1.5,'CapSize',0)
plot([0;1],log10([plsLoss(1);plsLoss(1)]),'k--','LineWidth',1.5)
% shadow the suboptimal entropy region
ah = gca;
% entRange = range(summH(:,1));
thd = 0.95;  % range of the entropy
TOL = 1e-3;
cs = spaps((0.05:0.05:1)',summH(:,1),TOL,1,3);
d1 = fnder(cs);
minSp = fnzeros(d1,[0.1,0.99]);
maxEntr = fnval(cs, minSp(1,1));
minEntr = summH(1,1);
refEntr = minEntr + (maxEntr - minEntr)*thd;
subOptRange = fnzeros(fncmb(cs,'-',refEntr),[0.1,0.99]);

% add two vertical lines  show the optimal entropy
yR = ah.YLim;
plot(ah,[subOptRange(1,1);subOptRange(1,1)],yR','--','LineWidth',1.5)
plot(ah,[subOptRange(1,2);subOptRange(1,2)],yR','--','LineWidth',1.5)
box on
set(gca,'YLim',[-0.65,0.15])
xlabel('$\rho_w$','interpreter','latex')
ylabel('$\log_{10}$(error)','interpreter','latex')


% add extral three lines to show sp = 0.1, 0.55,
plot(ah,[0.1;0.1],yR','--','LineWidth',1.5)
plot(ah,[0.6;0.6],yR','--','LineWidth',1.5)
plot(ah,[0.95;0.95],yR','--','LineWidth',1.5)
hold off


figNamePref = ['recon_tensor',num2str(N),'M',num2str(M),...
    'sp',num2str(spar),'_nSig',num2str(noise),'_',date];
saveas(gcf,[saveFolder,filesep,figNamePref,'.fig'])
print('-depsc',[saveFolder,filesep,figNamePref,'.eps'])

%% Figure 6, classification and entropy
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/decoding';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% classData = 'classify_svm_N100M10H500_sp3_G2_ns0.05_24-Sep-2018.mat'; % this is SVM
% classData = 'classify_LDA_N100M10diffPattern_HS500_sp3_group2_13-Sep-2018.mat';
classData = 'LDA_N100M10sp3_ns0.05_nType2_P200_H500_dS0.1_26-Dec-2018.mat';

entrData = 'entropyDistr_N100M10_sp3_sig2_ns0.05_26-Sep-2018.mat';
% load(fullfile(dFolder,classData))
load('/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding/LDA_N100M10_diffCluster/LDA_N100M10sp3_ns0.05_nType2_P100_H500_dS0.1_27-Dec-2018.mat')
load(fullfile(dFolder,entrData))

% parameters
N = 100;
M = 10;
sp = 3;
sig = 2;
noiseSig = 0.05;

% allSp = 0.05:0.05:1;
allSp = 0.1:0.05:1;

% selected data
ix = 4;   % corresponds to 100 patterns
% rawData = squeeze(summData.errorRate(:,:,ix));
rawData = allTestError;

fh = figure;
hold on
Bu = brewermap(11,'Blues');    % blues
errorbar(allSp',mean(rawData,2),std(rawData,0,2),...
        'o-','MarkerSize',10,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:),'LineWidth',1.5,'CapSize',0)
% lg = legend('N = 100,M=10,\sigma_n = 0.05, P = 100');
lg = legend('classification');
legend boxoff
set(lg, 'FontSize',16)

xlabel('$\rho_w$','interpreter','latex')
ylabel('classification error')

% add the second y axis to show the differential entropy
% yyaxis right
% ylabel('differential entropy','Color',Gr(8,:))
% errorbar(allSp(1:19)',summH(:,1),summH(:,2),'d-','MarkerSize',8,'Color',Gr(8,:),...
%     'MarkerFaceColor',Gr(8,:),'LineWidth',1,'CapSize',0)

ah = gca;
% y2h = ah.YAxis(2);
% y2h.Color = 'k';
box on


% shadow the suboptimal entropy region
entRange = range(summH(:,1));
thd = 0.95;  % range of the entropy
TOL = 1e-3;
cs = spaps((0.05:0.05:0.95)',summH(:,1),TOL,1,3);
% figure
% hold on
% fnplt(cs)
% scatter(allSp(1:19)',summH(:,1),'o')
% hold off
% first order derivative
d1 = fnder(cs);
minSp = fnzeros(d1,[0.1,0.99]);
maxEntr = fnval(cs, minSp(1,1));
minEntr = summH(1,1);
refEntr = minEntr + (maxEntr - minEntr)*thd;
subOptRange = fnzeros(fncmb(cs,'-',refEntr),[0.1,0.99]);

% add two vertical lines 
yR = ah.YLim;
plot(ah,[subOptRange(1,1);subOptRange(1,1)],yR','--','LineWidth',1.5)
plot(ah,[subOptRange(1,2);subOptRange(1,2)],yR','--','LineWidth',1.5)
hold off

figNamePref = ['class_LDA_N',num2str(N),'M',num2str(M),...
    'sp',num2str(sp),'_nSig',num2str(noiseSig),'_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% Figure 7 Both excitation and inhibition
% load the summary data
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

dFile = 'gcmi_inhiActi_diffBasal_diffN_QY_24-Jun-2018.mat';
load(fullfile(dFolder,dFile))

N = 50;  %in the main figure only plot one N
inx = find(nRecp ==N);

% we also need some simulation on when r0 = 0;
load([dFolder,filesep,'gcmi_ref',filesep,'int_N50_R10_S2_sig2_2018-07-17.mat'])
f0 = mean(-allfmin);
stdf0= std(-allfmin);

% ========================================================================
% plot how differential entorpy change with r_0
% ========================================================================
figure
% hold on
errorbar(allR0',meanFmin(:,inx),stdFmin(:,inx),'o-','MarkerSize',...
    12,'Color',Bu(10,:),'MarkerFaceColor',Bu(10,:),'LineWidth',2,'CapSize',0)
% errorbar(0,f0,stdf0,'d','MarkerFaceColor',myRd,'MarkerEdgeColor',myRd)
% hold off
xlabel('$r_0$','Interpreter','latex','FontSize',28)
ylabel('differential entropy','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear','XTick',0:0.25:1)
figNamePref = ['Figure6_ExciInhi_fmin_basal_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ========================================================================
% plot how inhibitory fraction changes with r0
% ========================================================================
figure
hold on
errorbar(allR0',allRatio(:,inx),stdRatio(:,inx),'o-','MarkerSize',...
    12,'Color',Bu(10,:),'MarkerFaceColor',Bu(10,:),'LineWidth',2,'CapSize',0)
plot([0;1],[0;1],'--','LineWidth',1.5,'Color',Gr(7,:))
hold off

box on
xlabel('$r_0$','Interpreter','latex','FontSize',28)
ylabel('inhibitory fraction','FontSize',28)
set(gca,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['Figure6_ExciInhi_inhiFrac_basal_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% ========================================================================
% plot the example histogram
% ========================================================================
% exampFile = 'int_N50_R10_S2_sig2_alp0.18_2018-06-22.mat';
% load([dFolder,filesep,'gcmi_inhi',filesep,exampFile])

selFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_inhi/N50M10S2sig2_1013';
exampFile = fullfile(selFolder,'gcmiInhi_N50_R10_S2_sig2_alp0.18_frac_2018-10-13.mat');
load(exampFile)

% essential parameters
nOdor = 50;
nRecp = 10;
sig = 2;
sp = 2;
r0 = 0.18;

inx = 1;
w = allMat(:,1);
we = log(w(allSign(:,1)>0));
wi = log(w(allSign(:,1)<0));

% histogram
figureSize = [0 0 5 4.5];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
[ha, pos] = tight_subplot(2,1,[0 0],[.2 .01],[.16 .03]);

axes(ha(1))
h1 = histfit(we,15);
xlim([-4,3.5])
ylim([0,100])
set(gca,'XTick',[],'YTick',0:50:100,'FontSize',24);

axes(ha(2))
h2 =  histfit(wi,15);
xlim([-4,3.5])
ylim([0,100])
xlabel('$\ln(w)$','Interpreter','latex')
set(gca,'XTick',-4:2:2,'YTick',0:50:100,'FontSize',24);

% set the color
h1(1).FaceColor = Bu(5,:);h2(1).FaceColor = RdBu(4,:);
h1(2).Color = Bu(10,:);h2(2).Color =RdBu(2,:);
h1(2).LineWidth = 3;h2(2).LineWidth = 3;

figNamePref = ['Figure6_ExciInhi_example_histo_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ========================================================================
% plot the  heatmap
% ========================================================================
allBlue = brewermap(128,'Blues');
allRd = brewermap(128,'Reds');

% w = reshape(w,[nRecp,nOdor]);
sign = reshape(allSign(:,1),[nRecp,nOdor]);

[inx_Y,inx_X] = find(sign >0);
[inx_Yinhi,inx_Xinhi] = find(sign <0);
% symbolExci = zeros(nRecp,nOdor);

%map value to color
LB = -2;
UB = 3;
we(we < LB) = LB;
we(we > UB) = UB;

LB_inhi = -3;
UB_inhi = 2;
wi(wi < LB_inhi) = LB_inhi;
wi(wi > UB_inhi) = UB_inhi;
% colorExci = zeros(length(we),3);
colorExci = allBlue(round((we - LB)/(UB-LB)*127)+1,:);
colorInhi = allRd(round((wi - LB_inhi)/(UB_inhi- LB_inhi)*127)+1,:);

maxMarker = 20;
allSizeExci = (we - LB)/(UB-LB)*(maxMarker - 1)+1;
allSizeInhi = (wi - LB_inhi)/(UB_inhi-LB_inhi)*(maxMarker - 1)+1;

figure
hold on
for i0 = 1:length(we)
    plot(inx_X(i0)-0.5,inx_Y(i0)-0.5,'o','MarkerSize',allSizeExci(i0),'MarkerFaceColor',...
        colorExci(i0,:),'MarkerEdgeColor',colorExci(i0,:),'LineWidth',0.5)
end

for i0 = 1:length(wi)
     plot(inx_Xinhi(i0)-0.5,inx_Yinhi(i0)-0.5,'o','MarkerSize',allSizeInhi(i0),'MarkerFaceColor',...
        colorInhi(i0,:),'MarkerEdgeColor',colorInhi(i0,:),'LineWidth',0.5)
end
hold off
set(gca,'XTick',10:10:50,'YTick',[1,5,10])
% grid on
% set(gca,'XTick',[])
% daspect([1 1 1])
figNamePref = ['Figure6_ExciInhi_examp_heatmap_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% Plot how optimal W parameters changes with N, M, n, but using paramterized W distribution
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

dName = 'gmci_distri_summData_29-Jul-2018.mat';
load(fullfile(dFolder,dName));

% ============================================================
% plot how target function changes with sp for differnt M
% ============================================================
% we only plot N = 100
figure
hold on
for i0 = 1:length(allM) - 1
    errorbar(allSp',allMeanFmin{1}(i0,:)',allStdFmin{1}(i0,:)','o-','MarkerSize',...
    12,'Color',Bu(1 + 2*i0,:),'MarkerFaceColor',Bu(1 + 2*i0,:),'LineWidth',2,'CapSize',0)
end
hold off
xlabel('$n$','Interpreter','latex')
ylabel('$H$','Interpreter','latex')
box on
figNamePref = ['Suppl_GcmiDistr_fmin_n_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% ============================================================
% plot how std changes with sp for differnt M
% ============================================================

% we only plot N = 100
figure
hold on
for i0 = 1:length(allM) - 1
    errorbar(allSp',allMeanSigW{1}(i0,:)',allStdSpW{1}(i0,:)','o-','MarkerSize',...
    12,'Color',Bu(1 + 2*i0,:),'MarkerFaceColor',Bu(1 + 2*i0,:),'LineWidth',2,'CapSize',0)
end
hold off
box on

xlabel('$n$','Interpreter','latex')
ylabel('$\sigma_{w}$','Interpreter','latex')

figNamePref = ['Suppl_GcmiDistr_sigW_n_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% ============================================================
% plot how sparsity of W changes with sp for differnt M
% ============================================================

% we only plot N = 100
figure
hold on
for i0 = 1:length(allM) - 1
    errorbar(allSp',allMeanSpW{1}(i0,:)',allStdSpW{1}(i0,:)','o-','MarkerSize',...
    12,'Color',Bu(1 + 2*i0,:),'MarkerFaceColor',Bu(1 + 2*i0,:),'LineWidth',2,'CapSize',0)
end
ld = legend('M=5','M=10','M=15','M=20','M=30','Location','northeast');
set(ld,'FontSize',16)
legend boxoff
hold off
xlabel('$n$','Interpreter','latex')
ylabel('$\rho_{w}$','Interpreter','latex')
box on
figNamePref = ['Suppl_GcmiDistr_rhoW_n_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% ============================================================
% plot how mean of W changes with sp for differnt M
% ============================================================

% we only plot N = 100
figure
hold on
for i0 = 1:length(allM) - 1
    errorbar(allSp',allMeanAveW{1}(i0,:)',allStdAveW{1}(i0,:)','o-','MarkerSize',...
    12,'Color',Bu(1 + 2*i0,:),'MarkerFaceColor',Bu(1 + 2*i0,:),'LineWidth',2,'CapSize',0)
end
hold off
box on
xlabel('$n$','Interpreter','latex')
ylabel('$\mu_{w}$','Interpreter','latex')

figNamePref = ['Suppl_GcmiDistr_muW_n_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% Suppl: scatter plot overlay with histogram to show the uniform distribution of response

% use my simulation 09/06/2018
% dFolder = '/Users/shan/Dropbox/olfactionProject/data/gcmi';
% dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_Ndp';

N = 80;
M = 13;
sp = 3;
sig = 2;

% need to set a break point
optMatrixCMA_v6(N,M,sp,sig);

% load the slected matrix as optimal W
% dFolder = '/Users/shan/Dropbox/olfactionProject/data/GcmiN100R30-10-03';
% dFolder = '/Users/shan/Dropbox/olfactionProject/data/fig3/GcmiNdp1018';
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_Ndp';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% fName = 'N100_R30_S2_sig2_2018-10-03.mat';
fName = 'Gcmi_N80_R13_S3_sig2_2018-10-21.mat';

load(fullfile(dFolder,fName))
w0 = allMat(:,1);

% now run the specific line in "optMatrixCMA_v6", but make sure set the
% break point in "gcmiCost"
[wmin,fmin,~,~,~,bestever] = cmaes('gcmiCost',w0,iniSig,opts);

%then plot the scatter pont and histogram
ix1 = 4;
ix2 = 10; 

figure('Renderer', 'Painters')
subplot(1,3,1)
scatter(resp(ix1,:),resp(ix2,:),'.')
xlabel('ORN1','FontSize',20)
ylabel('ORN2','FontSize',20)

subplot(1,3,2)
hold on
h1 = histogram(resp(ix1,:),'Normalization','pdf','FaceColor',lBu,'FaceAlpha',0.4,...
    'EdgeColor','none');
stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
    [0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',2)
xlabel('response of ORN1','FontSize',20)
ylabel('probability','FontSize',20)
hold off
box on

subplot(1,3,3)
hold on
h2 = histogram(resp(ix2,:),'Normalization','pdf','FaceColor',lBu,'FaceAlpha',0.4,...
    'EdgeColor','none');
stairs([h2.BinEdges(1),h2.BinEdges,h2.BinEdges(end)],...
    [0,h2.Values,h2.Values(end),0],'Color',dpBu,'LineWidth',2)
xlabel('response of ORN2','FontSize',20)
ylabel('probability','FontSize',20)
hold off
box on

prefix = ['suppl_scatter',date];
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])

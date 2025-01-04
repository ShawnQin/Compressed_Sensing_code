% reconfiguration of the figures, t

close all
clear

%% load the original fiure, default in 
figFolder = '/Users/shan/Dropbox/olfactionProject';
subFolder = 'poster figure';
figResPath = fullfile(figFolder,subFolder,'sigmaSens.fig');
expResFig = openfig(figResPath);
expResFig.WindowStyle = 'normal';
expResFig.DockControls = 'off';
expResFig.PaperPositionMode = 'auto';

%% set the parameter of graphics
figureSize = [0,0,3.4,2.5];
axisFontSize = 14;
tickWidth = 1;
labelFontSize = 14;
colorSet = [3,110,184;224,135,51;202,39,44;0,0,0]/256;


%% reset the figure
expResFig.PaperPositionMode = 'auto';
expResFig.Units = 'inches';
% expResFig.Position = figureSize;
expResFig.CurrentAxes.XLim = [0,800];
expResFig.CurrentAxes.XTick = 0:200:800;
set(expResFig.CurrentAxes,'LineWidth',tickWidth,'FontName','Helvetica',...
    'FontSize',axisFontSize)
xlabel('Time(a.u)','FontName','Helvetica','FontSize',labelFontSize)
ylabel('Residule','FontName','Helvetica','FontSize',labelFontSize)
%set(gca,'XLabel','Time','YLabel','Residule','FontSize',labelFontSize)
figNamefig = fullfile(figFolder,subFolder,'timeDetrExpfprPRE.fig');
figNamepdf = fullfile(figFolder,subFolder,'timeDetrExpfprPRE.eps');
saveas(expResFig,figNamefig)
expResFig.PaperUnits = 'inches';
%expFig.PaperPosition = figureSize;
%print(expResFig,'-painters','-depsc',figNamepdf)
expResFig.Position = figureSize;
print(expResFig,'-depsc',figNamepdf)
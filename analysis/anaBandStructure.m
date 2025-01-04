% This program test wheter w_1 and w_2 on the the ring of 
% band structure of optimal W matrix produce similar output (entropy).
% To put it in another way, if this has some degeneration

close all
clear


% load the parameters that determine the band structure
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
fileName = 'summData_0418_Regul.mat';
data = load(fullfile(dFolder,fileName));

%%%% define my personal color platte
re = [];  %red
gr = [];  %green
bl = [];  %blue
cy = [];  %cyan
or = [];  %orange
gy = [];  %grey
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
set(groot, 'defaultLineLineWidth',3);
set(groot, 'DefaultAxesColor', 'remove');
set(groot, 'DefaultAxesLineWidth', 1.5);
set(groot, 'DefaultFigureInvertHardcopy', 'on');
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.2);
set(groot, 'DefaultAxesFontSize', 24);
set(groot, 'DefaultLegendBox', 'off');
set(groot, 'DefaultLegendFontSize',20);

set(groot,'DefaultFigurePaperUnits', 'inches', ...
           'DefaultFigureUnits', 'inches', ...
           'DefaultFigurePaperPosition', [0, 0, 6.1, 4.6], ...
           'DefaultFigurePaperSize', [6.1, 4.6], ...
           'DefaultFigurePosition', [0.1, 0.1, 6, 4.5]);

greyColors = brewermap(11,'Greys');
blues = brewermap(11,'Blues');
%%
% set the parameters
global param trainData

param.nOdor = 2;                    % number of odors
param.numR = 100;                   % number of receptors
param.lmu= 0;  
param.lSig= 2.9;                    % std for lognorm
param.h = 1;                        % Hill coefficient of response function
param.regularize = false;           % if apply regularization on 
numSamp = 1e4;
    
trainData = mvnrnd([param.lmu, param.lmu],[param.lSig^2,0;0,param.lSig^2],numSamp);
%determine the p and c0 depending of sigma_c
index_sig = find(data.allSig == param.lSig);
index_R = find(data.R == param.numR);
p = data.p{index_R,index_sig}(1);
c0 = data.c0{index_R,index_sig}(1);

% Random select some points on the band curves and estimates the response entropy
% define random select w1 and calcualte w2
numPoints = 100;
w1 = -rand(numPoints,1)*3;
w2 = -abs(c0 - abs(w1).^p).^(1/p);

allH = zeros(numPoints,1);


% ecdf
figure
hold on
for i0 = 1:numPoints
    W = [w1(i0);w2(i0)]*param.lSig;
    resp = exp(trainData)*exp(W)./(1 + exp(trainData)*exp(W));
    [f, x] = ecdf(resp);
    plot(x,f,'Color',greyColors(6,:),'LineWidth',1.5)
    mult = 1;
    col = HShannon_KDP_initialization(mult);
    allH(i0) = HShannon_KDP_estimation(resp',col);
end
hold off
title(['$\sigma_c =',num2str(param.lSig),',|w_1|^{',num2str(p),'} + |w2|^{',...
    num2str(p),'} = ',num2str(c0),'$'],'Interpreter','latex')
xlabel('response')
ylabel('ecdf of response')
figPref = ['ecdf_w1_w2_onband',date];
% saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
% print('-depsc',[saveFolder,filesep,figPref,'.eps'])


% distribution of differential entropy
figure
histogram(allH,'Normalization','probability')
title(['$\sigma_c =',num2str(param.lSig),',|w_1|^{',num2str(p),'} + |w2|^{',...
    num2str(p),'} = ',num2str(c0),'$'],'Interpreter','latex')
xlabel('differential entropy')
ylabel('probability')
figPref = ['pdf_H_w1_w2_onband',date];
% saveas(gcf,[saveFolder,filesep,figPref,'.fig'])
% print('-depsc',[saveFolder,filesep,figPref,'.eps'])

figure
scatter(w1,w2)
title(['$\sigma_c =',num2str(param.lSig),',|w_1|^{',num2str(p),'} + |w2|^{',...
    num2str(p),'} = ',num2str(c0),'$'],'Interpreter','latex')
% lg = legend(['$\sigma_c =',num2str(param.lSig),',|w_1|^{',num2str(p),'} + |w2|^{',...
%     num2str(p),'} = ',num2str(c0),'$']);
% set(lg,'Interpreter','latex')
xlim([-3 0])
ylim([-3 0])
xlabel('$w_1/\sigma_c$','Interpreter','latex')
ylabel('$w_2/\sigma_c$','Interpreter','latex')

figure
histogram(resp,'Normalization','probability')
xlabel('response')
ylabel('probability')


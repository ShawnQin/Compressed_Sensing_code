% this script analyze the PN response data
% the first data set is from Badel et al, Neuron, 2016 using Ca2+ imaging
% method. The data contains PN response to both pure and mixture of odors
% we want to see the statistics of PN response, and compare it with OSN

% last revised on 07/10/2018

clear
clc

%%  load and prepare the data
dFolder = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/data';
file = 'Badel_Neuron2016.xls';

[NUM,TXT,RAW] = xlsread(fullfile(dFolder,file));

% the first  35 rows are for pure odorants, the rest are for mixture

outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

file = 'Badel_Neuron2016.xls';
[NUM,TXT,RAW] = xlsread(fullfile(dFolder,file));

inx0 = 34; % first 34 rows corresponds to pure odorant
% the first  35 rows are for pure odorants, the rest are for mixture
% number of strong repsonsive glomeruli (excitation or inhibition) for each
% odorant
strongGlo = [5,14,18,6,11,14,0,2,7,6,4,17,5,10,11,10,19,13,13,3,4,8,13,11,...
    12,5,1,4,2,4,3,5,8,2,13,9,2];

% number of odors one glomerulus response to (strong interaction)
strongOdor = [19,8,18,6,8,15,21,21,3,9,3,9,20,16,9,6,8,3,6,5,6,6,9,6,0,4,1,...
    5,8,0,8,4,3,13,3,3,2];

%% estimate the EC50 or Km
% here we assume a simple PN response as a function of odor concentration
% Hill function
n0 = 1;  % we can do better by fitting the dose response curve
Rmax = ceil(max(max(NUM)));  % the maximum response
temp = NUM(1:inx0,:)/Rmax;
allKm = 1./temp(temp>0.001)-1;

%% plot the figures
defaultGraphicsSetttings

% how the avearge response of each glomerulus
figure
bar(sort(mean(NUM/100,1),'descend'))
xlabel('glomeruli index')
ylabel('averge response ($\Delta F/F$)','Interpreter','latex')
figNamePref = ['avePNresp_BadelNeuron2016_hist_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% histogram of all PN response
figure
histogram(NUM(:)/100,'Normalization','probability')
xlabel('response ($\Delta F/F$)','Interpreter','latex')
ylabel('probability')
figNamePref = ['PNresp_BadelNeuron2016_hist_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])
<<<<<<< HEAD
=======

% histogram of all PN response for pure odorants
figure
temp1 = NUM(1:inx0,:);
histogram(temp1(:)/100,'Normalization','probability')
xlabel('response ($\Delta F/F$)','Interpreter','latex')
ylabel('probability')
figNamePref = ['PNrespPure_BadelNeuron2016_hist_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% histogram of all PN response for odor mixture with different
% concentration
figure
temp2 = NUM(inx0+1:end,:);
histogram(temp2(:)/100,'Normalization','probability')
xlabel('response ($\Delta F/F$)','Interpreter','latex')
ylabel('probability')
figNamePref = ['PNrespMixture_BadelNeuron2016_hist_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% histogram of strong glomeruli response
figure
histogram(strongGlo)
xlabel('number of glomeruli')
ylabel('Frequency')

% histogram of number of odorants that a glomerulus strongly response
figure
histogram(strongOdor)
xlabel('number of odor')
ylabel('Frequency')

%% effective Km of PN
figure
histogram(log(allKm))
xlim([-2,6])
xlabel('effective $\ln(K_m)$','Interpreter','latex')
ylabel('frequencey')
figNamePref = ['effecPN_Km_exci_BadelNeuron2016_hist_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% correlation of the PN response pattern
% only use the pure odorant data
temp = NUM(1:inx0,:);
C = corr(temp');
imagesc(C)

% cluster
Y = pdist(temp,'spearman'); 
% Z = linkage(Y,'single');
% T = cluster(Z,'Cutoff',4); 
% dendrogram(Z,'ColorThreshold',4)

Z = linkage(temp','single','spearman');
[H,T,outperm] = dendrogram(Z);

C1 =  corr(temp(:,outperm));
imagesc(C1)

temp = NUM(1:inx0,:);
C = corr(temp);
imagesc(C)

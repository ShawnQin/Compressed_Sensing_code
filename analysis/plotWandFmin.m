% this program plot the figures during algorithm test

close all
clear

%% plot the distribution of optimal W and fmin

% data source
folder = '/Users/shan/Dropbox/olfactionProject/data/gcmi/gcmi_corr';

N = 20;
M = 9;
sp = 3;
sig = 2;

% allFNames = dir(folder);
% files = {allFNames.name}';
% 
% sFlag = '(?<= *R100_)';
% 
% for i0 = 1:length(sFlag)
%     
% end

% fileName1= 'IG_N1_R200_S2_sig2_2018-03-03.mat';
% fileName2= 'IG_N2_R100_S2_sig2_2018-03-06.mat';
% fileName3= 'IG_N3_R100_S3_sig2_2018-03-06.mat';
% fileName4= 'IG_N4_R100_S4_sig2_2018-03-06.mat';
% fileName5= 'IG_N5_R100_S5_sig2_2018-03-05.mat';  % N = 5
% fileName1= 'IG_N20_R9_S2_sig2_2018-03-09_noReg.mat';  % N = 5
fileName1= 'gcmi_N20_R9_S3_sig2_2018-03-20.mat'; 
% fileName1= 'IG_N20_R5_S3_sig2_2018-03-13.mat';


tp1 = load([folder,filesep,fileName1]);
% tp2 = load([folder,filesep,fileName2]);
% tp3 = load([folder,filesep,fileName3]);
% tp4 = load([folder,filesep,fileName4]);
% tp5 = load([folder,filesep,fileName5]);

% plot the distribution of W, N = 1
figure
histogram(tp1.allMat,'Normalization','probability')
xlabel('$\ln(W)$','Interpreter','latex','FontSize',28)
ylabel('probability','FontSize',28)
xlim([-50,10])
set(gca,'LineWidth',1.5,'FontSize',24)
legend(['N=',num2str(N),',M=',num2str(M),',sp=',num2str(sp),...
    ',\sigma=',num2str(sig)],'Location','northwest')
legend boxoff

figFlag = ['gcmi_hist_',num2str(N),'_M',num2str(M),'_sp',...
    num2str(sp),'_sig',num2str(sig),'_corr_',date];

figName=['./figureData/',figFlag];
saveas(gcf,figName)


% plot the W distribution of N=2
figure
histogram(tp1.allMat(:),'Normalization','probability')
xlabel('$\ln(W)$','Interpreter','latex','FontSize',28)
ylabel('probability','FontSize',28)
xlim([-40,10])
set(gca,'LineWidth',1.5,'FontSize',24)
legend(['N=',num2str(N),',M=',num2str(M),',sp=',num2str(sp),...
    ',\sigma=',num2str(sig)],'Location','northwest')
legend boxoff

figFlag = ['gcmi_hist_zoom_N',num2str(N),'_M',num2str(M),'_sp',...
    num2str(sp),'_sig',num2str(sig),date];
figName=['../figureData/','cmpIGwithInteg_histW_N2_M100_sig2_0303.fig'];
saveas(gcf,figName)

% check the active part of W
temp = tp1.allMat;
histogram(temp(abs(temp) <=10)/sig,'Normalization','probability')
xlabel('$\ln(W)/\sigma_c$','Interpreter','latex','FontSize',28)
ylabel('probability','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
xlim([-3,3])

legend(['N=',num2str(N),',M=',num2str(M),',sp=',num2str(sp),...
    ',\sigma=',num2str(sig)],'Location','northwest')
legend boxoff

figFlag = ['gcmi_hist_zoom_N',num2str(N),'_M',num2str(M),'_sp',...
    num2str(sp),'_sig',num2str(sig),'_corr_',date];
figName=['./figureData/',figFlag];
saveas(gcf,figName)

% individual matrix
sltInx = 2;
temp2 = tp2.allMat(:,sltInx);
histogram(temp2(temp2>-8))
xlabel('$\ln(W)$','Interpreter','latex','FontSize',28)
ylabel('count','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
legend(['N=',num2str(N),',M=',num2str(M),',sp=',num2str(sp),...
    ',\sigma=',num2str(sig)],'Location','northwest')
legend boxoff

% plot histogram of N=5

% check the active part of W
temp = tp3.allMat;
% histogram(temp(abs(temp) <8))
histogram(temp(:),'Normalization','probability')
xlim([-50,10])

xlabel('$\ln(W)$','Interpreter','latex','FontSize',28)
ylabel('count','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
legend('N= 5,M = 100,\sigma = 2','Location','northwest')
legend boxoff

% individual matrix
sltInx = 5;
temp3 = tp3.allMat(:,sltInx);
histogram(temp3(temp3>-8))
% xlim([-50,10])
xlabel('$\ln(W)$','Interpreter','latex','FontSize',28)
ylabel('count','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
legend('N= 5,M = 100,\sigma = 2','Location','northwest')
legend boxoff


% distribution of target function
figure
histogram(tp1.allfmin,'Normalization','probability')
xlabel('$f_{min}$','Interpreter','latex','FontSize',28)
ylabel('probability','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
xlim([0,3])

legend(['N=',num2str(N),',M=',num2str(M),',sp=',num2str(sp),...
    ',\sigma=',num2str(sig)],'Location','northwest')
legend boxoff

%% Heatmap of the individual matrxi
sltInx = 3;
temp = tp1.allMat(:,sltInx);
imagesc(reshape(temp,[M,N]),[-6 8])

%% summary of different W and N
allW = cell(5,1);
allW{1} = tp1.allMat;
allW{2} = tp2.allMat;
allW{3} = tp3.allMat;
allW{4} = tp4.allMat;
allW{5} = tp5.allMat;

allFmin = cell(5,1);
allFmin{1} = tp1.allfmin;
allFmin{2} = tp2.allfmin;
allFmin{3} = tp3.allfmin;
allFmin{4} = tp4.allfmin;
allFmin{5} = tp5.allfmin;

allSparsity = zeros(4,40);
allActive = zeros(5,40);
thd = -8;  %threshold is set as 4*sig
for i0 = 1:4
    allSparsity(i0,:) = sum(allW{i0+1}>thd,1)/((i0+1)*100);
end
meanSp = mean(allSparsity,2);
errSp = std(allSparsity,0,2);


% active position and std
meanActive = zeros(5,40);
stdActive = zeros(5,40);
for i0 = 1:5
    for j0 =  1:40
        temp = allW{i0}(:,j0);
        meanActive(i0,j0) = mean(temp(temp>thd));
        stdActive(i0,j0) = std(temp(temp>thd));
    end
        
end
aveMeanActive = mean(meanActive,2);
stdMeanActive = std(meanActive,0,2);

aveStdActive = mean(stdActive,2);
stdStdActive = std(stdActive,0,2);

%plot the N-dependent sparsity
figure
errorbar((2:5)',meanSp,errSp,'o-','MarkerSize',12,'LineWidth',2)
xlabel('number of odor','FontSize',28)
ylabel('sparsity of W','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
legend('M=100,\sigma =2')
legend boxoff

%plot the N-dependent active elements and position
figure
errorbar((1:5)',aveMeanActive,stdMeanActive,'o-','MarkerSize',12,'LineWidth',2)
xlabel('number of odor','FontSize',28)
ylabel('$$\langle W_{act}\rangle$$','Interpreter','latex','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
legend('M=100,\sigma =2')
legend boxoff

%plot the N-dependent active elements and std
figure
errorbar((1:5)',aveStdActive,stdStdActive,'o-','MarkerSize',12,'LineWidth',2)
xlabel('number of odor','FontSize',28)
ylabel('std. of $$\langle W_{act}\rangle$$','Interpreter','latex','FontSize',28)
set(gca,'LineWidth',1.5,'FontSize',24)
legend('M=100,\sigma =2')
legend boxoff


%% check the distribution of N = 2
dPath = '/Users/shan/Dropbox/olfactionProject/data/twoByMLargeSig';
fileName = '2xM_sig5_2018-03-02.mat';
load(fullfile(dPath,fileName));

figure
rawMat = allMat{15};
histogram(rawMat(:))

temp = rawMat(:);
figure
histogram(temp(temp > -8),40)

%scatter plot of pair of interactions
figure
hold on
for i0 = 1:size(rawMat,2)
    temp = reshape(rawMat(:,i0),[200,2]);
    scatter(temp(:,1),temp(:,2))
end
hold off

% example matrix
figure
inx = 3;
% eMat = reshape(rawMat(:,inx),[2,200]);
eMat = reshape(rawMat(:,inx),[200,2]);
imagesc(eMat,[-6,3])

figure
% scatter(eMat(1,:),eMat(2,:))
scatter(eMat(:,1),eMat(:,2))
xlabel('w1','FontSize',24)
ylabel('w2','FontSize',24)
set(gca,'FontSize',20,'LineWidth',1.5)


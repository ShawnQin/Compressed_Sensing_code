% this program read the optimized interaction matrix 
% and decompose it into two parts, one is "sensitive" and the other one 
% is no-response

% OUTPUT:
% allMoments  the first and second moment of each Kij matrix
% allNonzero  percentage of non-zero elements of each Kij matrix
% mu, sig     average and std of fitted Gaussian distribution of Ps(K)
% ps          percentage of non-zero elements of overall Kij
      
close all
clear
clc

% global varibles, used when fitting the lognormal distributed wij
% dataFolder = '/home/shan/Documents/MATLAB/theoNeurosci/olfaction/spLogN50_R5_S0.2_indep_whole';
% dataFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/data/spLogN50_R5_S0.1_indep_whole';
% allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
% files = {allFile.name}';
% 
% %merge data
% numRecp = 5;
% numOdor = 50;
% %prepare the data
% allMat = zeros(numOdor*numRecp,length(files));
% allfmin = zeros(length(files),1);
% for i0 = 1:length(files)
%     temp = load(char(fullfile(dataFolder,filesep,files{i0})));
%     allMat(:,i0)  = temp.wmin;
%     allfmin(i0) = temp.fmin;
% end

load('../data/indepMultiSpRatioWholeWithBoundsLowerN.mat');
nSparse = size(allWmin,2);
nRepeat = size(allWmin,1);
nRecep = 9;
nOdor = 20;
sp = 0.06:0.02:0.14;
allWmin = permute(allWmin,[3 1 2]);

for s = 1:nSparse
    
allMat = allWmin(:,:,floor(sp(s)/0.02-2));

%% fit each trial and estimate mu, sigma and p0
% only use data that has a absolute value smaller than 10
allMoments = zeros(size(allMat,2),2);  %first and second moment
allNonzero = zeros(size(allMat,2),1);  % percentage of non-zero elements

figure(1)
for i0 = 1:size(allMat,2)
    fitData = allMat(:,i0);
    temp = fitData(abs(fitData) <=10);
    subplot(nSparse,nRepeat,(s-1)*nRepeat+i0)
    histogram(temp,20,'Normalization','probability')
    hold on

    [f,xi] = ksdensity(temp,'width',1.2);
    plot(xi,f,'r--','LineWidth',1)
    [pks,locs] = findpeaks(f,xi);
    [~, Inx] = sort(pks,'descend');
    allMoments(i0,1) = locs(Inx(1));  % the highest peak is what we want
    threshold = exp(-1/2)*max(pks);
    %xleft = xi(find(f >= threshold,1,'first'));
    xright = xi(find(f >= threshold,1,'last'));
    allMoments(i0,2) = xright - allMoments(i0,1);  %estimated sigma of Gaussian distribution
    
    % find the position where to set the threhold of sensitive and
    % non-response Kij
%     [f2, x2] = ksdensity(fitData,'width',1.5);
%     [pks2,loc2]= findpeaks(-f2,x2);
%     [~, Inx2] = sort(pks2,'ascend');
    thresh = allMoments(i0,1) - 2*allMoments(i0,2);
    allNonzero(i0) = sum(sum(fitData >= thresh))/length(fitData(:));   
end

%% fit the overall distribution
temp = allMat(abs(allMat) <= 10);
[f3,x3] = ksdensity(temp(:),'width',1);
figure(2)
subplot(1,5,s)
plot(x3,f3,'b--','LineWidth',2)
hold on
histogram(temp(:),30,'Normalization','probability')
hold off
[pks3,locs3] = findpeaks(f3,x3);
[~, Inx3] = sort(pks3,'descend');
mu(s) = locs3(Inx3(1));
threshold = exp(-1/2)*max(pks);
xleft = x3(find(f3 >= threshold,1,'first'));
xright = x3(find(f3 >= threshold,1,'last'));
sig(s) = xright - mu(s);  %estimated sigma of Gaussian distribution

% estimate the nonzero Kij, use all the data
[f4,x4] = ksdensity(allMat(:),'width',1.5);
thresh = mu(s) - 2*sig(s);

figure(3)
subplot(1,5,s)
%% 
%% 
h = histogram(allMat(:),30,'Normalization','probability');
h.FaceColor = [1 1 1];
h.EdgeColor = 'k';
h.LineWidth = 1;
hold on
plot(x4,f4,'LineWidth',2)
hold off
legend('simulation','kernel fit')
legend boxoff
xlabel('ln(Kij)','FontSize',24)
ylabel('probability','FontSize',24)
set(gca,'XLim',[-25,20],'YLim',[0, 0.12],'LineWidth',1,'FontSize',20)

% [pks4,loc4]= findpeaks(-f4,x4);
% [~, Inx4] = sort(pks4,'ascend');
ps(s) = sum(sum(allMat >= thresh))/length(allMat(:));  % nonzeros percentage

end

nOdor = 20;
nRecep = 9;
nRepeat = 10;
allWmin = permute(allWmin,[4 1 2 3]);

%structure for discretized Kij
for t = 1:3
    for k = 1:3
        
    pattern = [];
    sp = 0.05*(t+1); 
    wmin = reshape(allWmin(:,:,t,k),nRecep,nOdor*nRepeat);
    wmindis = wmin;
%     thr = mu(floor(sp/0.02)-2)-2*sig(floor(sp/0.02)-2);
    thr = -10;
    wmindis(wmindis <= thr) = 0;
    wminnonzero = wmin(wmin > thr);
    part = prctile(wminnonzero,[1/3*100,2/3*100]); %set all Kij values smaller than -10 to zero
    wmindis(wmin>thr&wmin<part(1)) = 1;
    wmindis(wmin>part(1)&wmin<part(2)) = 2;
    wmindis(wmin>part(2)) = 3;
    wmindis = sort(wmindis,'descend');
    base = 4.^(0:8);
    %ps1 = ps(floor(sp/0.02)-2);
    ps1 = length(find(wmin<=-10))/1800;
    w = base*wmindis;
    stat = tabulate(w);
    num = stat(:,1);
    p = [];
    for i = 1:length(num)
        n = zeros(1,4);
        a = num(i);
        for j = 1:nRecep
        n1 = mod(a,4);
        a = fix(a/4);
        n(n1+1) = n(n1+1) + 1;
        end
        m(i,:) = n;
        p(i) = prod([nchoosek(nRecep,n(1)),nchoosek(nRecep-n(1),n(2)),nchoosek(nRecep-n(1)-n(2),n(3)),1].*[1-ps1,ps1/3,ps1/3,ps1/3].^n);
    end
    per = [stat(:,3)./100,p'];
    num = num(per(:,1)>0);
    per = per(per(:,1)>0,:);
    decimal = num(abs(per(:,2).^2-per(:,1).^2)>0.01);
    for i = 1:length(decimal)
         pat = zeros(1,4);
         dec = decimal(i);
    for j = 1:nRecep
        pat1 = mod(dec,4);
        dec = fix(dec/4);
        pat(pat1+1) = pat(pat1+1)+1;
    end
    pattern(:,i) = pat';
    end
    patternstat{t,k} = pattern; %obtain the patterns that are significant
%     figure(4)
% subplot(1,5,t)
% plot(per(:,1),'LineWidth',2)
% hold on
% plot(per(:,2),'LineWidth',2)
% xlabel('different responses')
% ylabel('probability')
% legend('simulation','theoretical')
figure(5)
subplot(3,3,3*(t-1)+k)
plot(per(:,1),per(:,2),'.','MarkerSize',25);
xlim([0,max(max(per))+0.005]);
ylim([0,max(max(per))+0.005]);
xlabel('simulation');
ylabel('random theoretical');
set(gca,'FontSize',14)
    end
end
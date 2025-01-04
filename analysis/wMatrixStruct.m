% this program using Guangwei's the 'GetBestOrder.m' to show if the
% optimizaed matrix is competable with his experiment

% 

close all
clear

%% choose the data to analyze
% dataFolder = '/Users/shan/Dropbox/olfactionProject/data/gcmi/gcmi_struct';
dataFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_struct';
fName = '50x20_indep_gcmi_sp5.mat';
load(fullfile(dataFolder,fName))
nOdor = 50;
nRecp = 20;
% w0 = reshape(allMat(:,1),[nRecp,nOdor]);
w0 = reshape(wmin,[nRecp,nOdor]);


figure
imagesc(w0,[-6,4])
xlabel('odorants')
ylabel('receptors')

% the correlation matrix after transformed
Z = corr(exp(w0));
figure
imagesc(triu(Z,1),[0,1])
set(gca,'FontSize',16)
xlabel('odorants','FontSize',20)
ylabel('odorants','FontSize',20)
%reference correlation matrix
D0 = allD0{1};
figure
imagesc(D0,[-1,1])
xlabel('odorants')
ylabel('odorants')
% rearrange hte order
[odorOrder, ORNOrder] = GetBestOrder(exp(w0));

%show the final rearranged matrix
w1 = w0(odorOrder,:);
w2 = w1(:,ORNOrder);
imagesc(w2',[-8,4])
xlabel('receptors')
ylabel('odorants')

% w3 = w0(:,randperm(100,20))';
w3 = w0(:,1:20)';

[odorOrder, ORNOrder] = GetBestOrder(exp(w3));
w4 = w3(odorOrder,:);
w5 = w4(:,ORNOrder);
imagesc(w5,[-6,4])
xlabel('odorants')
ylabel('odorants')

% first, find the strongest odorants corresponding to 
[X,inx] = sort(w0,2,'descend');
wselect = w0(:,inx(:,1));
[odorOrder, ORNOrder] = GetBestOrder(exp(wselect));
w6 = wselect(odorOrder,:);
w7 = w6(:,ORNOrder);
imagesc(w7,[-6,4])
set(gca,'FontSize',16)
xlabel('odorants','FontSize',20)
ylabel('odorants','FontSize',20)

%% calculate the correlation matrix of Guangwei's data

% load the data
fileName = '/Users/shan/Documents/MATLAB/theoNeurosci/larval_olfaction-master';
dataName = 'EC50Km.mat';
load(fullfile(fileName,dataName))

% transform the K matrix into correlation matrix
c0 = 4;  % the reference concentration
inx_nonzero = ec50Map > 0;
newW = zeros(size(ec50Map));
newW(inx_nonzero) = exp(ec50Map(inx_nonzero)/c0)./(1+exp(ec50Map(inx_nonzero)/c0));

corrMRecp = corr(newW);
corrMOdor = corr(newW'); %revised on 6/11/2018

save('./figureData/EC50CorrOdorMatrix.mat','corrMOdor')
save('./figureData/EC50CorrRecpMatrix.mat','corrMRecp')
% plot the correlation among odorants
figure
imagesc(corrMOdor)
colorbar
set(gca,'FontSize',20)
xlabel('odorants','FontSize',24)
ylabel('odorants','FontSize',24)


% plot correlation among receptors
figure
imagesc(corrMRecp)
colorbar
set(gca,'FontSize',20)
xlabel('receptors','FontSize',24)
ylabel('receptors','FontSize',24)

tm = corrMRecp;
tm(tm < 0.3) = nan;
imagesc(tm)
colorbar
set(gca,'FontSize',20)
xlabel('receptors','FontSize',24)
ylabel('receptors','FontSize',24)



tm = corrMRecp;
tm(tm > -0.3) = nan;
imagesc(tm)
colorbar
set(gca,'FontSize',20)
xlabel('receptors','FontSize',24)
ylabel('receptors','FontSize',24)


alp = 2;
X = rand(1e4,1);
xmin  = 1;
Y = xmin*(1-X).^(1/(-alp + 1));
[f1,x1] = ecdf(Y);
plot(log10(x1),log10(1-f1))

%% imposing correlation based on experiment data

% first, load the reference correlation matrix
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
corrFile = 'EC50CorrRecpMatrix.mat';
load(fullfile(dFolder,corrFile));
inxD = logical(triu(ones(size(corrMRecp)),1));

N = 18;  % number of odor
R = 18;  % number of receptor


%load the data with different strength of regularization
datafiles = dir([dFolder,filesep,'gcmi_struct',filesep,'*.mat']);
files = {datafiles.name}';

% define the variables to store the summarized data
allMat = zeros(N*R,length(files));
allfmin = zeros(length(files),1);
allDist = zeros(length(files),1);
allSp = zeros(length(files),1);    %input sparsity
allLbd = zeros(length(files),1);   %lambda of regularization

for i0  = 1:length(files)
    s1 = '(?<= *N)[\d.]+(?=_)';
    s2 = '(?<= *_R)[\d]+(?=_)';
    s3 = '(?<= *_S)[\d]+(?=_)';
%     s4 = '^\w{3,4}(?=N)';
    s5 = '(?<= *lambda)[\d.]+(?=\.mat)';
    numOdor = str2num(char(regexp(files{i0},s1,'match')));
    numRecp = str2num(char(regexp(files{i0},s2,'match')));
    allSp(i0) = str2num(char(regexp(files{i0},s3,'match')));
%     tp = char(regexp(files{i0},s4,'match'));
    allLbd(i0) = str2num(char(regexp(files{i0},s5,'match')));
    
    temp = load(char(fullfile(dFolder,filesep,'gcmi_struct',filesep,files{i0})));
    allMat(:,i0) = temp.wmin;
    allfmin(i0) = temp.fmin;
    
    w = reshape(temp.wmin,[numRecp,numOdor]);
    D = corr(exp(w')./(1+exp(w')));
    allDist(i0) = norm(D(inxD) - corrMRecp(inxD));
    
end

% order and rearrange the data
uniqSp = sort(uniq(allSp));
uniqLbd = sort(uniq(allLbd));

% best performance
bestInx = find(allDist == min(allDist));
w = reshape(allMat(:,bestInx),[R,N]);

figure
imagesc(corr(w))
[odorOrder, ORNOrder] = GetBestOrder(exp(w'));

% plot the reordered matrix
% w1 = w(odorOrder,:);
w2 = w(odorOrder,:);
figure
imagesc(corr(w2))

sum(corr(w) - corr(w2))

imagesc(w2,[-6,4])
set(gca,'FontSize',16)
xlabel('odorants','FontSize',20)
ylabel('odorants','FontSize',20)

figure
for i0 = 1:24
    subplot(4,6,i0)
    temp = reshape(allMat(:,i0),[18,18]);
    imagesc(temp,[-6,6])
end

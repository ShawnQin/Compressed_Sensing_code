% This program test how the response pattern of ORNs change with the
% concentrations. Several papers have suggested that odor identiy can be
% extracted if the relative pattern across different ORNs doesn't change
% very much

clear
clc

% color setting
RdBu = brewermap(11,'RdBu');   % red and blue
Bu = brewermap(11,'Blues');    % blues
Gr = brewermap(11,'Greys');    % greys

%% Excitation-only situtation
% load one of the optimal sensitivity matrix
dFolder = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/code/data/GcmiN100R30-10-03';
fName = 'N100_R30_S2_sig2_2018-10-03.mat';
load(fullfile(dFolder,fName));

% some basic parameters
N = 100;
M = 30;
sp = 2;
sig =2;

% select one matrix and show its sparsity
ix0 = 1;   % slect one matrix, default the first
w = reshape(allMat(:,ix0),[M,N]);  %the example optimal sensitivity matrix

% define the response function, input
respFun = @(X,w) exp(w)*X./(1 + exp(w)*X);

% genereate sample odors
inx = randi(N);  % random select one odorant
conRange = 10.^(-3:0.5:3);  % concentration range
X0 = zeros(N,length(conRange));
X0(inx,:) = conRange;

% response across different concentrations
resp = respFun(X0,w);

% plot the average response versus the concentration
figure
plot(conRange,mean(resp,1),'o-','MarkerSize',12,'LineWidth',2)
set(gca,'XScale','log','YScale','log','XTick',10.^(-3:1:3))
set(gca,'FontSize',24,'LineWidth',1.5)
xlabel('$\log_{10}(c)$','Interpreter','latex','FontSize',28)
ylabel('$\langle r \rangle$','Interpreter','latex','FontSize',28)


%% mixture with different concentration
% since the optimal W is determined undeter n = 2, we need assign the odor
% mixture to have the same sparsity

% generate odor mixtures
conRange = 10.^(-2:0.5:2);
vecC = exp(sig*randn(sp,1))*conRange;  % random concentration
inx = randi(N,[2,1]);    % random assign the position
C = zeros(N,size(vecC,2));
C(inx,:) = vecC;

resp = respFun(C,w);

% plot the average response versus the concentration
figure
plot(conRange,mean(resp,1),'o-','MarkerSize',12,'LineWidth',2)
set(gca,'XScale','log','YScale','log','XTick',10.^(-3:1:3))
set(gca,'FontSize',24,'LineWidth',1.5)
xlabel('$\log_{10}(c)$','Interpreter','latex','FontSize',28)
ylabel('$\langle r \rangle$','Interpreter','latex','FontSize',28)

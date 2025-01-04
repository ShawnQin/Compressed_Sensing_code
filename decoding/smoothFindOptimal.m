% this program smooth the performance v.s the sparsity of W and locate the
% position of optimal performance
% last revised on 09/19/2018

close all
clear

clear

% load the file 
dFoler = './';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
fName = 'classify_LDA_N100M10diffPattern_HS500_sp3_group2_13-Sep-2018.mat';
load(fullfile(dFoler,fName))

% basic parameters of this data set
N = 100;
M = 10;
sp = 3;
sig = 2;
H = 500;   % hidden layer size

[NS,NR,NP] = size(summData.errorRate);  % number of different patterns
optParam = zeros(NP,2);

for i0 = 1:NP
%     X = kron(ones(NR,1),summData.allSp');
    figure
    hold on
    Y = mean(summData.errorRate(:,:,i0),2);
    plot(summData.allSp',Y,'o')
    
    % cubic spline fit
    if i0 == 1
        TOL = 2e-4;
    else
        TOL = 1e-5;
    end
    cs = spaps(summData.allSp',Y,TOL,1,3);
    fnplt(cs)
    hold off
    
    % first order derivative
    d1 = fnder(cs);
    minSp = fnzeros(d1,[0.1,0.99]);
    optParam(i0,1) = minSp(1,1);
    optParam(i0,2) = fnval(cs,optParam(i0,1));
end

selInx = 2:NP;

% setting of graphics
defaultGraphicsSetttings
Bu = brewermap(11,'Blues');    % blues
seqColor = brewermap(21,'YlGnBu');


figure
plot(summData.allPattern(selInx)',optParam(selInx,1),'o-','LineWidth',2,'MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:))
ylim([0.2,0.6])
set(gca,'XScale','log')
xlabel('number of clusters')
ylabel('$\rho_w^*$','Interpreter','latex')
figNamePref = ['classif_N',num2str(N),'M',num2str(M),'sp',num2str(sp),'_Hidden',...
    num2str(H),'optSpW_diffCluster_',date];
saveas(gcf,[sFolder,filesep,figNamePref,'.fig'])
print('-depsc',[sFolder,filesep,figNamePref,'.eps'])


figure
plot(summData.allPattern(selInx)',optParam(selInx,2),'o-','LineWidth',2,'MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:))
ylim([0.001,0.3])
xlabel('number of clusters')
ylabel('best performance')
set(gca,'XScale','log')
figNamePref = ['classif_N',num2str(N),'M',num2str(M),'sp',num2str(sp),'_Hidden',...
    num2str(H),'bestPerf_diffCluster_',date];
saveas(gcf,[sFolder,filesep,figNamePref,'.fig'])
print('-depsc',[sFolder,filesep,figNamePref,'.eps'])



figure
hold on
for i0 = 1:length(selInx)
    plot(summData.allSp',mean(summData.errorRate(:,:,selInx(i0)),2),'Color',seqColor(3 + 3*i0,:))
%     [ix, iy] = sort(meanError(:,i0));
%     bestPerf(i0,1) = summData.allSp(iy(1));
%     bestPerf(i0,2) = ix(1);
end
lg = legend('cluster = 20','cluster = 50','cluster = 100','cluster = 200','cluster = 500',...
    'cluster = 1000','Location','eastoutside');
set(lg,'FontSize',16)
% legend box off
hold off
xlabel('$\rho_w$','Interpreter','latex')
ylabel('error')
figNamePref = ['classif_N',num2str(N),'M',num2str(M),'sp',num2str(sp),'_Hidden',...
    num2str(H),'diffCluster_',date];
saveas(gcf,[sFolder,filesep,figNamePref,'.fig'])
print('-depsc',[sFolder,filesep,figNamePref,'.eps'])


%% Find the smallest error for tensorflow reconstruction with respect number of hidden layer
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding/N100M20ns05diffL_1014';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
allFile = dir(fullfile(dFolder,filesep,'*.mat'));
files = {allFile.name}';

% parameters
N = 100;
M = 20;
spar = 2;
sig = 2;
noise = 0.05;
H = 500;

allSp = 0.05:0.05:1;
allL = 1:1:5;

repeat = 40;
summError = zeros(length(allSp),length(allL),repeat);
summMeanError = zeros(length(allSp),length(allL));
summStdError = zeros(length(allSp),length(allL));
optParam =  zeros(length(allL),2);

for i0 = 1:length(files)
%     st = '(?<= *_sp)[\d*]+(?=_)';
    ly = '(?<= *_L)[\d]+(?=.mat)';
    ss = '(?<= *_sp)[\d.]+(?=_)';
    layer = str2num(char(regexp(files{i0},ly,'match')));
    sp = str2num(char(regexp(files{i0},ss,'match')));
    ix = find(round(100*allSp) == round(100*sp));
   
    load(char(fullfile(dFolder,filesep,files{i0})));
    summError(ix,layer,:) = testLost;
    summMeanError(ix,layer) = mean(testLost);
    summStdError(ix,layer) = std(testLost);
    
end

% fit and find the best performance, by fitting a smooth curve
ix = 1:1:18;
for i0 = 1:length(allL)
    X = kron(ones(repeat,1),(1-allSp)');
    
    figure
    hold on
%     Y = squeeze(summError(:,i0,:));
%     plot(X', log(Y(:)),'.')
    
    plot((1-allSp(ix))',log(summMeanError(ix,i0)),'o')
    
    % cubic spline fit
    if i0 == 5
        TOL = 1e-2;
    else
        TOL = 2e-3;
    end
%     cs = spaps(X', log(Y(:)),TOL,1,3);
    cs = spaps((1-allSp(ix))', log(summMeanError(ix,i0)),TOL,1,3);
    fnplt(cs)
    hold off
    
    d1 = fnder(cs);
    minSp = fnzeros(d1,[0.1,0.99]);
    optParam(i0,1) = minSp(1,1);
    optParam(i0,2) = exp(fnval(cs,optParam(i0,1)));

end

% save data
prefix = ['tsRecon_summ_N',num2str(N),'M',num2str(M),'sp',num2str(spar),'ns',...
    num2str(noise),'H', num2str(H),'diffL_',date,'.mat'];
save(fullfile(sFolder,prefix),'summMeanError','summStdError','M','N','spar',...
    'noise','H','allSp','allL','sig','optParam')


% plot the figures

% for i0 = 1:length(allL)
% %     X = kron(ones(NR,1),summData.allSp');
%     figure
%     hold on
%     plot((1- allSp)', summMeanError(:,i0),'o')
%     
%     % cubic spline fit
%     if i0 == 1
%         TOL = 2e-4;
%     else
%         TOL = 1e-5;
%     end
%     cs = spaps((1- allSp)',summMeanError(:,i0),TOL,1,3);
%     fnplt(cs)
%     hold off
%     
%     % first order derivative
%     d1 = fnder(cs);
%     minSp = fnzeros(d1,[0.1,0.99]);
%     optParam(i0,1) = minSp(1,1);
%     optParam(i0,2) = fnval(cs,optParam(i0,1));
% end

% setting of graphics
defaultGraphicsSetttings
Bu = brewermap(11,'Blues');    % blues
seqColor = brewermap(21,'YlGnBu');

figure
plot(allL', optParam(:,2),'o-','MarkerSize',12,'MarkerFaceColor',Bu(9,:),...
    'MarkerEdgeColor',Bu(9,:),'LineWidth', 2)
ylim([0.2 1])
xlabel('layers')
ylabel('minimum error')
figNamePref = ['tensorRecons_N',num2str(N),'M',num2str(M),'sp',num2str(spar),'_Hidden',...
    num2str(H),'diffL',date];
saveas(gcf,[sFolder,filesep,figNamePref,'.fig'])
print('-depsc',[sFolder,filesep,figNamePref,'.eps'])


figure
plot(allL', optParam(:,1),'o-','MarkerSize',12,'MarkerFaceColor',Bu(9,:),...
    'MarkerEdgeColor',Bu(9,:),'LineWidth', 2)
ylim([0, 1])
xlabel('layers')
ylabel('$\rho_w^*$', 'Interpreter','latex')
figNamePref = ['tensorRecons_N',num2str(N),'M',num2str(M),'sp',num2str(spar),'_Hidden',...
    num2str(H),'diffL_rhow',date];
saveas(gcf,[sFolder,filesep,figNamePref,'.fig'])
print('-depsc',[sFolder,filesep,figNamePref,'.eps'])

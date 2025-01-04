clc
clear
load('indepMultiSpRatioCon_no_noise.mat')
global param 

param.nOdor = 50;           %number of total odors
param.nRecep = 5;           %number of receptors

param.dType = 'lognorm';  %odor concentratio distribution
param.method = 'wholeMat';   %optimization method, whole matrix or distribution, "wholeMat/distr"
param.nSamp = 1e4;        %number of odor mixture

param.lMu = 0;   % mean for lognorm
param.lSig= 2;    % std for lognorm
param.gMu = 5;   %mean for gaussian data
param.gSig = 1;  %std for gaussian data

param.sparsity = true;   %sparse input
param.corr = false;       %correlation across odor concentration
param.regul = false;     %if use further reuglarization in cost function
param.noise = true; %if add noise
param.noiseSig = 0.1;
%param.noiseSig = noise(k);
param.eig = [3,3,2,1];%specify the first 4 component in data,column vectors
param.spRatio = 0.04;     %sparse ratio, percentage of non-zero elements
eigVal = specifyEig(param.nOdor,param.eig);
corrCoefMat = randCorrCoef('buildin',eigVal);
testData = genTrainData(param,corrCoefMat);

figure('color','white')
for i = 1:9
subplot(3,3,i)
histogram(allWmin{i,3,1},20);ylim([0,150]);
xlabel('ln(weight)')
ylabel('probability distribution')
title(strcat('weight statistics with sparsity   ',num2str(i*0.02+0.02)))
end
saveas(gcf,'fcon_50by5_indep_sparsity0.003to0.009_gaussiannoise0.1.fig')

w = reshape(allWmin{5,3,1},5,50);
col1 = HShannon_KDP_initialization(1);
col2 = HShannon_vME_initialization(1);
resp = exp(w) * testData./(1+ exp(w) * testData);
if param.noise
resp = resp + normrnd(0,param.noiseSig,size(resp));
resp(resp>1)=1;
resp(resp<0)=0;
end
h1 = HShannon_KDP_estimation(resp,col1);
h2 = HShannon_vME_estimation(resp,col2);
disp(strcat('h1 - h2 =',32,num2str(h1-h2)));

figure('color','white');
for i = 1:4
subplot(2,2,i)
scatter(resp(i,:),resp(i+1,:),'o')
xlabel('first receptor');
ylabel('second receptor');
title('scatter plot of two receptor')
end
saveas(gcf,'test for uniformity f.fig')



figure('color','white')
for i = 1:5 
for j = 1:5
subplot(5,5,i+(j-1)*5)
histogram(allWmin{2*j-1,i,1},20)
xlabel('ln(weight)')
ylabel('count');ylim([0,size(allWmin{2*j-1,i,1},1)/2]);
title(strcat('odor number-',num2str(i*10+20),32,'sparsity-',num2str(0.04*j)))
end
end
saveas(gcf,'fcon_indep_sparsity0.003to0.009_nodor30to70_gaussiannoise0.1.fig')

figure('color','white')
for i = 1:3
for j = 1:5
subplot(3,5,j+(i-1)*5)
histogram(allWmin{j*2-1,3,i},20)
xlabel('ln(weight)')
ylabel('count');
ylim([0,150])
title(strcat('noise sigma-',num2str(i*0.1),32,'sparsity-',num2str(0.04*j)))
end
end
saveas(gcf,'fixedcon_50by5_indep_sparsity0.003to0.009_gaussiannoise0.1to0.3.fig')

figure('color','white')
color = ['r','y','g','b','m'];
linetype = ['*','o','.'];
for i = 1:5
   for j = 1:3
        plot(0.04:0.02:0.2,allFmin(:,i,j),strcat(linetype(j),'-',color(i)),'LineWidth',1)
        hold on
   end
end
xlim([0.04,0.2]);
title('allFmin');
box off;
set(gca,'FontSize',16,'LineWidth',1);
saveas(gcf,'allFmin_con_f_smallnoise')

figure('color','white')
for i = 1:5
    h = subplot(5,1,i);
    wmin = reshape(allWmin{2*i-1,3},5,50);
    wmin(wmin<=-10) = -inf;
    wmin(wmin>=10) = inf;
    wmin = [wmin;zeros(1,50)];
    wmin = [wmin,zeros(6,1)];
    pcolor(wmin);
    caxis([-10,2]);
    title(strcat('odor = 50, sparsity =',32,num2str(0.04*i)));
    PN = get(h,'pos');
    PN = PN + [0.05 0 0 0];
    set(h,'pos',PN);    
end
colormap jet
bar = colorbar;
set(bar,'Position',[0.05 0.05 0.05 0.9])
saveas(gcf,'f_wmin_property_50by5.fig')

allMat = permute(allWmin,[]);
for i = 1:4
    for j = 1:5
  subplot(i,j,5*(i-1)+j)
  histogram(allMat(:,:,i,j));
  xlabel('ln(weight)');
  ylabel('count');
  ylim([0,200]);
    end
end


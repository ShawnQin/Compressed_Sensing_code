% This program make animation of how the optimal matrix evolve with time
% We are trying to make three movies: the first one is the heatmap of marix
% the second one is the hitogram of all elements, and the third one is the
% histogram of active elements
% 

close all
clear

%%
%load the file
% dataFolder = '/home/shan/Documents/MATLAB/theoNeurosci/olfaction/debug_N40_R16_S3_sig2_IT5000';
dataFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_animation/spLogN20_R9_S3_sig2_24';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

filenameprefix = 'outcmaes';
filenameextension = '.dat';
objectvarname = 'xrecentbest';
raw = load([dataFolder,filesep,filenameprefix objectvarname filenameextension]);
allWmin = raw(:,6:end);

allf= load([dataFolder,filesep,filenameprefix 'fit' filenameextension]);  %all fitness
fitness = allf(:,6);  %target function values for each iteration
iter = raw(:,1);     %iteration number

%% set up the first frame
% figureSize = [4 4 8 4];
set(groot,'default')

figure('Color','white')
[ha, pos] = tight_subplot(2,1,[.1,.1],[.2,.1],[.1,.05]);
% set(gcf,'Units','pixels','Position',figureSize,...
% 'PaperPositionMode','auto');
% set the parameters related to
nRecp = 9;
nOdor = 20;
sp = 3;
sig = 2;

t = reshape(allWmin(1,:),[nRecp,nOdor]);

% subplot(2,1,1)
axes(ha(1))
fh1 = imagesc(t,[-6,6]);
colorbar
ht = title(sprintf('Iteration: %d',1));


% subplot(2,1,2)
axes(ha(2))
fh2 = semilogy(iter(1),fitness(1),'b','LineWidth',3);
axis([0 iter(end) 1 1e8])
xlabel('Iteration','FontSize',16)
ylabel('Loss function','FontSize',16)
set(gca,'LineWidth',1,'FontSize',14)

% get figure size
set(gcf,'Units','pixels')
pos =get(gcf,'Position');
width = pos(3);
height = pos(4);
%% loop through all matrix

% preallocate data
% T = size(allWmin,1);
totIter = 2e3;  %we only need the first 1e3 iterations
allInx = [1:1:20,22:2:100,110:350];  %select part of the matrix
% allInx = 1:5:size(allWmin,1); % there is total 1e4 iteration
T = length(allInx);

% CAUSION: the factor 2 happens only for retina display !!!
mov = zeros(2*height,2*width,1,T,'uint8');

for i0 = 1:T
    % update graphics data. This is more efficient than recreatig plots
    fh1.CData = reshape(allWmin(allInx(i0),:),[nRecp,nOdor]);
    ht.String = ['Iteration: ',num2str(iter(allInx(i0)))];
    
    fh2.XData = iter(1:allInx(i0));
    fh2.YData = fitness(1:allInx(i0));
    % Get frame as an image
    f = getframe(gcf);
    
    if i0 == 1
        [mov(:,:,1,i0),map] = rgb2ind(f.cdata,256,'nodither');
    else
        mov(:,:,1,i0) = rgb2ind(f.cdata,map,'nodither');
    end
    
end

%Create animated GIF
movName = fullfile(saveFolder,['W_N',num2str(nOdor),'M',num2str(nRecp),'sp',...
    num2str(sp),'sig',num2str(sig),'_Reg.gif']);
imwrite(mov,map,movName,'DelayTime',0,'LoopCount',inf)


%% sparsity change of the matrix
spW = zeros(T,1);
for i0 = 1:T
    spW(i0) = sum(allWmin(allInx(i0),:)>-4*2)/nRecp/nOdor;
end
defaultGraphicsSetttings
figure
plot(iter(allInx),spW)
set(gca,'XScale','log','XTick',[1,10,1e2,1e3,1e4])
xlabel('iterations')
ylabel('sparsity of W')
set(gca,'FontSize',24,'LineWidth',1.5)
figNamePref = ['sparsityW_iteration_N',num2str(nOdor),'_M',num2str(nRecp),...
    '_sp',num2str(sp),'_sig',num2str(sig),'_',date];
saveas(gcf,[saveFolder,filesep,figNamePref,'.fig'])
print('-depsc',[saveFolder,filesep,figNamePref,'.eps'])

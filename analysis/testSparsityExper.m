% This program try to demonstrate that the undetected sensitivity matrix is
% essentially zero even though they could be very small. And this cannot be
% simply explained by the experiment measurement accuracy

%% color setting
defaultGraphicsSetttings
RdBu = brewermap(11,'RdBu');   % red and blue
Bu = brewermap(11,'Blues');    % blues
Gr = brewermap(11,'Greys');    % blues

lBu = [96,166,223]/255; %light blue

%% load the data
file3 = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/flyLarva_Guangwei_LogEC50_New.xlsx';
saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

[NUM3,TEXT,~]=xlsread(file3);
allW = NUM3(abs(NUM3) > 0);

newM = -NUM3;

%% Fit the data with a skewed Gaussian distribution
gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
skewStand = @(x,alpha) 2*gaussian(x).*normcdf(alpha*x);
skewPdf =  @(x,xi,omi,alp) 2/omi*gaussian((x-xi)/omi).*normcdf(alp*(x-xi)/omi);

dVec = newM(~isnan(newM));
% ecdf
[Y,X] = ecdf(dVec);

fun = @(x,xdata) normcdf((xdata-x(1))/x(2)) - 2*myOwenT((xdata-x(1))/x(2),x(3));
x0 = [mean(dVec),std(dVec),3];
lb = [0,1,0];
ub = [10,10,10];
% least square fit
optParam = lsqcurvefit(fun,x0,X(2:end),Y(2:end),lb,ub);

xi = optParam(1); omi = optParam(2); alp = optParam(3);


% ================================================
% fit histogram and fit with Gaussian distribution
% ================================================
figure
figureSize = [0 0 5 4.5];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
hold on
set(gcf,'renderer','Painters')
hh = histogram(abs(-allW),15,'Normalization','pdf',...
    'EdgeColor','none','FaceColor',lBu,'FaceAlpha',0.4);

X = 0:0.05:10;  % The range
Y2 = skewPdf(X,xi,omi,alp);
plot(X,Y2,'LineWidth',3,'Color',RdBu(2,:))
box on
ylim([0,0.5])
hold off

lg = legend('experiment','Skewed Gaussian fit');
set(lg,'FontSize',16)
% legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
set(gca,'Layer', 'top')


%% Assume a detection function and derive the "real" pdf of sensitivity
% we assume a sigmoidal function of detection threshold
thd = log(10^2);
n = 5;
detectFun = @(x) 1./(1+exp(-n*(x-thd)));
% X = -3:0.1:4;
% figure
% plot(X,detectFun(X))

% putative probability density function of sensitivity
hiddenPdf = @(x) skewPdf(x,xi,omi,alp)./detectFun(x);
X = 0:0.05:10;
figure
plot(X,hiddenPdf(X))
plot(X,skewPdf(x,xi,omi,alp))


%% 
n = 1.42;  % Hill coefficient
r = 0.05;
logW = @(r,n) 4 + 1/n*log10(r./(1-r));
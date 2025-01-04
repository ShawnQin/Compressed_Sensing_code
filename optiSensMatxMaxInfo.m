function [wmin,fmin] = optiSensMatxMaxInfo(numOdor,numRecp,sig,varargin)
%========================================================================
%This program try to find the optimal sensitvitiy matrix for a olfactory
%system using gradient descent method. The target function is the
%mutual information between input and output 
%The current program only works for overcomplete situation!

%writen by Shanshan Qin, Qianyi Li,Tang Lab

%VERSION 0.3  
%     add weight decay (L2 norm) to the target function to eliminate too
%     large W. Update of W is done in the logscale
%last revised 1/28/2018
%========================================================================
%{
p = inputParser;
p.addRequired('numOdor',@isnumeric);
p.addRequired('numRecp',@isnumeric);
% p.addRequired('sparsity',@isnumeric);
p.addRequired('sig',@isnumeric);

p.addParameter('numSamp',1e3,@isnumeric);  % maximum iteration times
p.addParameter('maxIter',2e3,@isnumeric);  % maximum iteration times
p.addParameter('learnRate',1,@isnumeric);  % learning rate
p.addParameter('repeats',1,@isnumeric);   % repeat times

parse(p,numOdor,numRecp,sig,varargin{:})

numOdor = p.Results.numOdor;
numRecp = p.Results.numRecp;
numSamp = p.Results.numSamp;
inputSig = p.Results.sig;
numRepeat = p.Results.repeats;
eta = p.Results.learnRate;
maxIter = p.Results.maxIter;
%}
%%

numSamp = 1e4;
inputSig = sig;
numRepeat = 20;
eta = 20;
maxIter = 2e3;

%% ============================================================
% global allInput

% when run the program on cluster, add the path to this file 
addpath(genpath('/lustre1/tangc_pkuhpc/ssqin/olfactoryCoding/optiInterMatrix'));
% warning('off','all')  % turn off the warning message
%% define basic parameters of the olfactory system'
% numRecp = 20;     %number of receptors
% numOdor = 3;      %number of odorants

meanCon = 0;        %average concentration in log scale
% sigma = 2;        %define the width of Gaussian distribution

% eta = 5;          %learning rate
% batchSize = 1e3;  %batch size
% maxIter = 1e4;    %maximum iteration times
     
Threshold = 1e-5;   %if target function doesn't change
% numSamp = 1e3;    %total smapled odor mixture

numWorker = 10;    %number of workers
%matlabpool close force local
%c = parcluster();
matlabpool(numWorker);     % this only works for matlab13 or older version
%parpool(c);

%% define some anoymous functions
resFunMM = @(x) x./(1+x); %the M-M response function
resFunMMd1 = @(x) (1+x).^(-2); % first order derivative
resFunMMd2 = @(x) -2*(1+x).^(-3); %second order derivative

% resFunLg = @(x) 1./(1+exp(-x)); %logistic(sigmoid) response function
% resFunLgd1 = @(y) y.*(1-y);    %first order derivative
% resFunLgd2 = @(y) y.*(1-y).*(1-2*y); %second order

%% generate a sample of odor space
% MU = meanCon*ones(1,numOdor);
% SIGMA = sigma*ones(1,numOdor);
% corrCoefMatrx = normrnd(0,0.5,[numOdor,numOdor]);
% corrCoefMatrx(corrCoefMatrx>1) = 1;
% corrCoefMatrx(corrCoefMatrx<-1) = -1;
% corrCoefMatrx = 1/2*(corrCoefMatrx'+corrCoefMatrx);
% corrCoefMatrx(logical(eye(size(corrCoefMatrx))))= 1; %set the diagonal elements to be 1
% 
% sigMatx = SIGMA'*ones(1,numOdor);
% covMatx = corrCoefMatrx.*sigMatx.*sigMatx';
% allInput = exp(mvnrnd(MU,covMatx,numSamp));

% generate identical independent samples
allInput = cell(numWorker,1);
for i0 = 1:numWorker
    allInput{i0} = exp(normrnd(meanCon,inputSig,[round(numSamp/numWorker),numOdor]));
end
%% specify the system
% aveChg = zeros(numRecp,numOdor);
% costFun = zeros(numOdor,numOdor);

% numRepeat = 2;          %repeats number, start from different initial values
totalCost = cell(numRepeat,1);   %record all the cost function
allW = zeros(numRecp*numOdor,numRepeat);  %record all the final W matrix
fmin = zeros(numRepeat,1);
wmin = zeros(numOdor*numRecp,numRepeat);
lastVal = 0;     %used to store the last value of last time target function

for j0 =1:numRepeat
    allWeight = [];
    weight = min(exp(normrnd(0,3*inputSig,[numRecp,numOdor])),exp(4*inputSig));  %initalize weight matrix
    iterTime = 0;
    changeVal = Inf; %
    
    while (iterTime < maxIter && changeVal > Threshold)
%         miniBatch = allInput(randperm(numSamp,batchSize),:);
        [weight,cost] = updateFun(weight,eta,numRecp,numOdor,allInput) ;
        totalCost{j0,1} = [totalCost{j0,1},cost];
        allWeight = [allWeight,weight(:)];
        iterTime = iterTime + 1;
    
        changeVal = abs(cost - lastVal);  %estimate the increasement of target function
        lastVal = cost;
    
        if mod(iterTime,100) ==0
            disp(['iteration ',num2str(iterTime), ',target function:', num2str(cost)])
        end
    end
    allW(:,j0) = weight(:);
    fmin(j0) = totalCost{j0,1}(end);
    wmin(:,j0) = weight(:);
end

% save the data
save('infoMax.mat','fmin','allW')
matlabpool close

% shutdown parpool
%p = gcp;
%delete(p)

% plot the evolution of target function
% plot(1:length(totalCost),totalCost,'LineWidth',2)
% xlabel('Iterations','FontSize',28)
% ylabel('differential entropy','FontSize',28)
% set(gca,'FontSize',24,'LineWidth',1)

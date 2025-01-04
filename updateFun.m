function [newWeight,costFun] = updateFun(weight,eta,numRecp,numOdor,allInput,varargin)
% this function implemente the learning precedure
% weight    weight matrix before update
% sample    a batch sample used to update
% eta       learning rate
% ====================================
% global allInput
%% parameter about the learing rate
alpha = 1e-1;      %weight decay constant

%% define some anoymous functions
resFun = @(x) x./(1+x); %the M-M response function
resFund1 = @(y) 1./(1+y).^2; % first order derivative
resFund2 = @(y) -2*1./(1+y).^3; %second order derivative

%% 

% aveChg = zeros(numRecp,numOdor);
% costFun = 0;
allCostFun = zeros(length(allInput),1);


% newWeight = weight;
N = size(allInput{1},1);
% inputData = allInput(randperm(N,1e3),:);

% parfor loop
allChg = zeros(numRecp,numOdor,length(allInput));
parfor j0 = 1:length(allInput)
    partData = allInput{j0};
for i0 = 1:N

%updation rule
vec = partData(i0,:);
totAct = weight*vec'; %total activation
fd1 = resFund1(totAct); %first order derivative
fd2 = resFund2(totAct); %second order derviative

%G matrix
gMx = diag(fd1);

%Kai matrix
kMx = gMx*weight + 1e-10*randn(size(gMx,1),size(weight,2)); %revised on
% kMx = gMx*weight; %revised on

%gamma matrix
%gmMx = pinv(kMx'*kMx)*kMx'*gMx;  %the computation of matrix inverse might be time consumming
% gmMx = (kMx'*kMx)\kMx'*gMx;  %the computation of matrix inverse might be time consumming
%gmMx = invChol_mex(kMx'*kMx)*kMx'*gMx;  %the computation of matrix inverse might be time consumming


% gamma vector
% gm = diag(kMx*gmMx).*fd2./fd1.^3;


%average over
% aveChg = aveChg + (gmMx'+gMx'*gm*vec)/size(allInput,1);
% costFun = costFun + 1/2*log(det(kMx'*kMx))/size(allInput,1);

% gmMx = pinv(kMx'*kMx)*kMx';
% gm = diag(kMx*gmMx).*fd2./fd1;

gmMx = pinv(kMx'*kMx);
gm = diag(kMx*gmMx*weight').*fd2;
allChg(:,:,j0) = allChg(:,:,j0) + (gMx*kMx*gmMx'+ gm*vec)/N;

% aveChg = aveChg + (gMx*gmMx'+ gm*vec)/size(allInput,1);
% costFun = costFun + 1/2*log(det(kMx'*kMx + 1e-20))/size(allInput,1);
% allChg(:,:,j0) = allChg(:,:,j0) + (gMx*gmMx'+ gm*vec)/N;

allCostFun(j0) = allCostFun(j0) + 1/2*log(det(kMx'*kMx + 1e-20))/N;




end
end
aveChg = mean(allChg,3);
costFun = mean(allCostFun);
%update weight matrix
% newWeight = weight + eta*aveChg;
% newWeight = exp(log(weight) + eta*weight.*aveChg);
% newWeight = exp(log(weight) + eta*weight.*aveChg) - eta*alpha*weight; %with weight decay
newWeight = max(exp(log(abs(weight)) + eta*abs(weight).*aveChg) - eta*alpha.*abs(weight).*heaviside(weight-exp(4*2)),0); %with weight decay

% newWeight = exp(log(weight) + ((1-alpha)*eta0 + alpha*eta)*weight.*aveChg);

end
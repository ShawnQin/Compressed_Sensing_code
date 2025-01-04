function MI = MI_QusiMonteCarlo(w)
% this function approximate mutual infomation between input and output with
% Gaussian noise, Using Quasi-Monte Carlo integration of I_G
% The quasi-random sequence can be generated through various methods
%
% w          the MxN weight matrix
% MI    an approximation of mutual information
%
% last revised on 02/28/2018

global param trainData


% get the relavent parameters from param
N = param.nOdor;                         % number of odor
sig_input = param.lSig;                  % std of input
sig_noise = param.noiseSig;              % std of noise
weight = reshape(w,[param.nRecep,N]);    % reshaped weight matrix


% precision matrix of input, used in I_G
Sm = 1/sig_input^2*eye(N);

%define Gaussian density function
mu = 0;
% gauss = @(x) 1/sqrt(2*pi)/sig_input*exp(-(x-mu).^2/sig_input^2);
% normpdf(x,0,sig)

% calculate the average value of I_G
costFun = 0;
Nsamp = param.nSamp;

% % generating quasi-random sequence,'sobolset' or 'lhsdesign'
% qrnd = quasiRandSeq(N,Nsamp,'lhsdesign');
% trainData = (qrnd - 0.5)*8*sig_input;

allNormDensity = normpdf(trainData,0,sig_input);
normIntegr = prod(allNormDensity*10*sig_input,1);
for i0 = 1:Nsamp

% vec = trainData(:,i0);
vec = exp(trainData(:,i0));  %revised on 03/04/2018


% this will improve the performance twice if the loop number is big
denomenator = (1+exp(weight)*vec).^2;
J = exp(weight).*vec'./denomenator;
fd1 = J'*J;

% fd1 = respFund1(weight,exp(vec)); %first order derivative

%G matrix
G = fd1 + Sm*sig_noise^2;   %debug on 02/28/2018
% G = fd1;   %debug on 02/28/2018

% costFun = costFun - 1/2*log(det(G))*prod(normpdf(vec,0,sig_input)*10*sig_input);
costFun = costFun - 1/2*log(det(G))*normIntegr(i0);

% costFun = costFun - 1/2*det(G);


if ~isreal(costFun)
    disp('complex vale of target function!')
end
% disp(costFun)
end

% averge value and adding the entropy of input to get the exact value of MI
% costFun = costFun/Nsamp + N/2*(log(2*pi) + 1);
% MI = costFun - N/2*log(2*pi*exp(1)*sig_input^2);

% % MI = costFun/Nsamp*10*sig_input;
MI = costFun/Nsamp;


end
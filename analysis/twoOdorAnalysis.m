% one odor, many receptors
% the seache    d wij value is already normalized by sigma, the width of input

function twoOdorAnalysis(sigma,varargin)
% sigma   an integer,the standard devivation of concentration distribution

% when run on cluster
addpath(genpath('/lustre1/tangc_pkuhpc/ssqin/olfactoryCoding/optiInterMatrix'));

N = 2; %number of odors
nRec = [1:1:20,30,40,50,100,200];
%nRec = [20 50];18-04-1

nsample = 1;   %for each sigma sample 50 times

allWmin = cell(length(nRec),1);
allFmin = nan(length(nRec),nsample);

for i0 = 1:length(nRec)
    allWmin{i0} = zeros(nRec(i0),nsample);
    for j0 = 1:nsample
        [wmin,fmin] = optMatrixCMA_v2_edited(N,nRec(i0),2,sigma);
        allWmin{i0}(:,j0) = wmin/sigma;   %normalized to sigma
        allFmin(i0,j0) = -fmin + log(2*pi*exp(1)*sigma*sigma);  % add input entropy
    end
end

save('oneByM.mat','allWmin','allFmin')
end

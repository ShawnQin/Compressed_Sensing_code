% this program compare the searched best interaction matrix using
% integration method
close all
clear

allSig = 0.5:0.5:10;
nsample = 20;   %for each sigma sample 50 times

allWmin = cell(length(allSig),1);
allFmin = nan(length(allSig),nsample);

for i0 = 1:length(allSig)
    allWmin{i0} = zeros(4,nsample);
    for j0 = 1:nsample
        [wmin,fmin] = optMatrixCMA_v2_edited(2,2,2,allSig(i0));
        allWmin{i0}(:,j0) = wmin/allSig(i0);   %normalized to sigma
        allFmin(i0,j0) = -fmin + log(2*pi*exp(1)*allSig(i0));  % add input entropy
    end
end

save('twobywo.mat','allWmin','allFmin')

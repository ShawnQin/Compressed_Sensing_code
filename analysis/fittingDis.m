%% a naive approach to fit distribution
allWmin = permute(allWmin,[1 3 2]);
for i = 1:4 
allMat = reshape(allWmin(:,:,i),1,6500);
[f,xi] = ksdensity(allMat,'width',1.2);
[pks,locs] = findpeaks(f,xi);[~, Inx] = sort(pks,'descend');
loc = locs(Inx(1));
threshold = exp(-1/2)*max(pks);
xright = xi(find(f >= threshold,1,'last'));
sigma = xright - loc;
pk(i) = loc;
sigmas(i) = sigma;xrights(i) = xright;
end


%% fitting the distribution via optimization
global param col
param.allMat = reshape(allWmin(:,:,2),1,size(allWmin,1)*size(allWmin,2));
col = DKL_kNN_k_initialization(1);
lambda = 2;
p0 = 0.3;
aver = 0;
sigma = 2;
pos = 1;
opts.PopSize = 30;
x0 = [lambda,p0,aver,sigma,pos];
opts.MaxIter = 1e3; %maximum number of iteration
opts.StopFitness = 0;
iniSig = [10,0.5,10,10,10]';
opts.LBounds = [0,0,-inf,0,-inf]';
opts.UBounds = [inf,1,inf,inf,inf]';
[wmin,fmin,~,~,~,bestever] = cmaes('costFunDis',x0,iniSig,opts);


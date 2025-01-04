function H = entropyCostDist(x)
% this function return differential entropy 
% interaction matrix is sampled from x, typcially the mean and std

% define global variables
global trainData  param
persistent  col mult
mult = 1;
col = HShannon_KDP_initialization(mult);


% repeats when estimate the differential entropy
N = 10;
Hlist = zeros(N,1);
for j0 =  1:N
    w = zeros(param.nRecep,param.nOdor);
    nonZeroInx = randperm(param.nRecep*param.nOdor,round(param.nRecep*param.nOdor*x(3)));
    w(nonZeroInx) = exp(normrnd(x(1),abs(x(2)),[length(nonZeroInx),1]));
    resp = w * trainData./(1+ w* trainData) + param.noiseSig*randn(param.nRecep,param.nSamp);


    MI = nonparanormal_info(resp');
    H0 = 0;
    for i0 = 1:param.nRecep
        H0 = H0 + HShannon_KDP_estimation(resp(i0,:),col);
    end
    Hlist(j0) = -H0 + MI;
end
H = mean(Hlist);
    %% regularization of the Weight
    if param.regularize
        lambda = 1;
        C = lambda*max((max(w(:))-5*param.lSig),0)*sum(sum(w(w>0).^2));
        H = H + C;
    end
end

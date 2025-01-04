function h = entropyCostMean(w)
% this is modified from entropyCost.m
% here a contraint on the average firing rate is imposed 
% This is used for optMatrixCMA_v4.m
% define persistent varible to save time
%
    global param X trainData
    %persistent trainData col mult
    mult = 1;
    col = HShannon_KDP_initialization(mult);

    w = reshape(w,[param.nRecep,param.nOdor]);
%% directly estimate the entropy                                       
    derive_resp = @(x) x.*(1-x);
    resp = exp(w) * trainData./(1+ exp(w) * trainData);
    if param.noise
       resp = resp + param.noiseSig.*normrnd(0,sqrt(derive_resp(resp))); %add noise
    end
    resp(resp>=1) = 1 - X(1:length(find(resp>=1)));  %avoid response to be equal to 1
    resp(resp<=0) = X(1:length(find(resp<=0))); 
        
    h = -HShannon_KDP_estimation(resp,col);

%%  define effective dimension as regularization

if param.regularize
    % two fold regularization, first we don't want the elements be small
    % we also want the average spiking rate be constant
    lambda = 1; %constraint on value of elements
    gm2 = 1;   %contraint coefficient on average rate
    C1 = lambda*max(((max(max(w)))-5*param.lSig),0)*sum(sum(w(w>0).^2));
    
    aveResp = mean(resp,2);
    C2 = gm2*norm(aveResp - ones(param.nRecep,1)*param.aveRate);
    
    h = h + C1 + C2;
end

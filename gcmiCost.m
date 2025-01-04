function H = gcmiCost(w)
% define persistent varible to save time
    global param trainData %trunctatePd
    
    mult = 1;
    col = HShannon_KDP_initialization(mult);
    
    w = reshape(w,[param.nRecep,param.nOdor]);
%     resp = exp(w)*trainData./(1+exp(w)*trainData) + random(trunctatePd,[param.nRecep,size(trainData,2)]);
    resp = exp(w)*trainData./(1+exp(w)*trainData) + param.noiseSig*randn(param.nRecep,size(trainData,2));
    
%     resp = exp(w)*trainData./(1+exp(w)*trainData);

    %% estimate the MI using gcmi
%     gcnomY = copnorm(resp');
%     MI = -mi_gg(gcnomX,gcnomY);
    MI = nonparanormal_info(resp');
%     MI = nonparanormal_Kendall_info(resp');
    %% using KDP to estimate the entropy of each dimension
    H0 = 0;
    for i0 = 1:param.nRecep
        H0 = H0 + HShannon_KDP_estimation(resp(i0,:),col);
    end
    H = -H0 + MI;   %without noise
%     H = -H0 + MI - param.nRecep/2*log(2*pi*exp(1)*param.noiseSig^2);
    
    %% regularization of the Weight
    if param.regularize
        lambda = 1;
%         C = lambda*max((max(w(:))-5*param.lSig),0)*sum(sum(w(w>0).^2));
        C = lambda*max((max(w(:))-5*param.lSig),0);
        H = H + C;
    end
end

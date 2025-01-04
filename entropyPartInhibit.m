function H = entropyPartInhibit(r)
% this fucntion return the differential entropy when partial receptors are
% considered to be inhibitory
% define persistent varible to save time
% r      the response of inhibitory neurons

    global param trainData 
    mult = 1;
    col = HShannon_KDP_initialization(mult);
    r_0 = param.r_0;
    alpha = (1-r_0)./r_0;
    
    % substitute 0 and 1 response to avoid singularity
    r(r==0) = 1e-5*abs(normrnd(0,1,size(r(r==0))));
    r(r==1) = 1-1e-5*abs(normrnd(0,1,size(r(r==1))));
    r = reshape(r,[param.nRecep,param.nOdor]);
    
    if param.frac == 1 %all inhibitory receptors
        Sign = 2*(r>r_0)-1;
        w = ((1./r-1)/alpha).^(-Sign)-1;

        % set the inhibitory and excitatory matrix
        w_in = zeros([param.nRecep,param.nOdor]);
        w_ac = w_in;
        w_in(Sign==-1) = w(Sign==-1);
        w_ac(Sign==1) = w(Sign==1);
        resp = (1+alpha*(1+w_in*trainData)./(1+w_ac*trainData)).^(-1) ...
            + param.noiseSig*normrnd(0,1,param.nRecep,param.nSamp);
    elseif param.frac == 0  % all excitatory
        w = 1./(1./r - 1);   
        resp = w*trainData./(1 + w*trainData) + param.noiseSig*normrnd(0,1,param.nRecep,param.nSamp);
        
    else
        % get the response of inhibitory receptors and positive-long receptors
        r1 = r(1:ceil(param.nRecep*param.frac),:);
        r2 = r(ceil(param.nRecep*param.frac)+1 : end,:);
        
        Sign = 2*(r1>r_0)-1;
%         w = ((1./r1-1)/alpha).^(-Sign)-1;
        w = exp(-Sign.*(log(1-r1) - log(alpha.*r1)))-1;  %modified by Qianyi
    
        % only positive response matrix
%         wp = 1./(1./r2 - 1);
        wp = exp(log(r2)-log(1-r2));  % modified by Qianyi

        % set the inhibitory and excitatory matrix
        w_in = zeros([ceil(param.nRecep*param.frac),param.nOdor]);
        w_ac = w_in;
        w_in(Sign==-1) = w(Sign==-1);
        w_ac(Sign==1) = w(Sign==1);
        resp = [(1+alpha.*(1+w_in*trainData)./(1+w_ac*trainData)).^(-1); ...
            + wp*trainData./(1+wp*trainData)] + param.noiseSig*normrnd(0,1,param.nRecep,param.nSamp);
    end
    %% estimate the MI using gcmi
    MI = nonparanormal_info(resp');
%     MI = nonparanormal_Kendall_info(resp');

    %% using KDP to estimate the entropy of each dimension
    H0 = 0;
    for i0 = 1:param.nRecep
        H0 = H0 + HShannon_KDP_estimation(resp(i0,:),col);
    end
    H = -H0 + MI;

    %% regularization of the Weight

    if param.regularize
        lambda = 1;
        ERR = 1e-4;
%         C = lambda*max(((max(max(log(w))))-5*param.lSig),0)*sum(sum(w(w>0).^2));
        C = lambda*(sum(sum(r1 > 1-ERR)) + sum(sum(r1< ERR)) + sum(sum(r2>1-ERR)));
        H = H + C;
    end
end
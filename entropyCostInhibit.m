function H = entropyCostInhibit(r)
% define persistent varible to save time
    global param trainData 
    mult = 1;
    col = HShannon_KDP_initialization(mult);
    r_0 = param.r_0;
    alpha = (1-r_0)/r_0;
    r = reshape(r,[param.nRecep,param.nOdor]);
    r(r==0) = 1e-5*abs(normrnd(0,1,size(r(r==0))));
    r(r==1) = 1-1e-5*abs(normrnd(0,1,size(r(r==1))));
    Sign = 2*(r>r_0)-1;
%     w = ((1./r-1)/alpha).^(-Sign)-1;
    w = exp(-Sign.*(log(1-r) - log(alpha*r)))-1;

    % set the inhibitory and excitatory matrix
    w_in = zeros([param.nRecep,param.nOdor]);
    w_ac = w_in;
    w_in(Sign==-1) = w(Sign==-1);
    w_ac(Sign==1) = w(Sign==1);
    resp = (1+alpha*(1+w_in*trainData)./(1+w_ac*trainData)).^(-1) ...
        + param.noiseSig*normrnd(0,1,param.nRecep,param.nSamp);
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
        C = lambda*max(((max(max(w)))-5*param.lSig),0)*sum(sum(w(w>0).^2));
        H = H + C;
    end
end
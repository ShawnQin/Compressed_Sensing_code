function h = entropyCostPure(w)
% define persistent varible to save time
% the cost function is just entropy, without any constriant and
% regularization
    global param trainData
    %persistent trainData col mult
    mult = 1;
    col = HShannon_KDP_initialization(mult);
    w = reshape(w,[param.nRecep,param.nOdor]);
%% directly estimate the entropy
    resp = exp(w) * trainData./(1+ exp(w) * trainData);
    h = -HShannon_KDP_estimation(resp,col);

%% discrete the response state to 0, 1, 2,3
   
%     bits = 4;
%     nr = param.nRecep;
%     base = bits.^(0:(nr-1));
%     r = base*floor(bits*exp(w) * trainData./(1+ exp(w) * trainData));  %on exponential scale
%     prob = tabulate(r);
%     nonzeop = prob(prob(:,3) > 0,3)/100;
%     h = sum(nonzeop.*log2(nonzeop));

end

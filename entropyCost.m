function h = entropyCost(w)
% define persistent varible to save time
    global param X trainData trunctatePd 
    %persistent trainData col mult
    mult = 1;
    if strcmp(param.target,'mi')
        col = IShannon_DKL_initialization(mult);
    else 
        col = HShannon_KDP_initialization(mult);
%         col = HShannon_kNN_k_initialization(mult);
%         col = HShannon_expF_initialization(mult);
    end
    w = reshape(w,[param.nRecep,param.nOdor]);
%% directly estimate the entropy
    if param.function
        resp = 1/2*(erf(log(exp(w)*trainData)/2)+1);
        resp(resp==1) = 1 - X(1:length(find(resp==1)));  %avoid response to be equal to 1
        resp(resp==0) = X(1:length(find(resp==0)));
    else                                       
%       derive_resp = @(x) x.*(1-x);
      resp = exp(w) * trainData./(1+ exp(w) * trainData);
      
%       if param.noise
%         resp = resp + param.noiseSig.*normrnd(0,sqrt(derive_resp(resp))); %add state dependent noise
%       end  
%      resp(resp>=1) = 1 - X(1:length(find(resp>=1)));  %avoid response to be equal to 1
%      
%      if any(resp<=0)
%         resp(resp<=0) = X(1:length(find(resp<=0))); 
%      end
    
      % add trucate Gaussian noise
     if param.noise           
%         pd = makedist('Normal','mu',0,'sigma',param.noiseSig);
%         trunctatePd = truncate(pd,0,0.2);
        resp = resp + random(trunctatePd,[param.nRecep,size(trainData,2)]);
     end
     
     if strcmp(param.target,'mi')
        h = -IShannon_DKL_estimation([resp;trainData],[param.nRecep,param.nOdor],col);
     else
        h = -HShannon_KDP_estimation(resp,col);
%         h = -HShannon_kNN_k_estimation(resp,col);
%         h = -HShannon_expF_estimation(resp,col);
     end
    end
%     h = -HShannon_vME_estimation(resp,col);
%     h = 0;
%     for i0 = 1:size(resp,1)
%         h = h  -HShannon_KDP_estimation(resp(i0,:),col);
%     end
%     
%     for i0 = 1:size(resp,1)
%         h = h - HShannon_vME_estimation(resp(i0,:),col2);
%     end
    %h =  - HShannon_kNN_k_estimation(exp(w) * trainData./(1+ exp(w) * trainData),col2);
%% discrete the response state to 0, 1, 2,3
   
%     bits = 4;
%     nr = param.nRecep;
%     base = bits.^(0:(nr-1));
%     r = base*floor(bits*exp(w) * trainData./(1+ exp(w) * trainData));  %on exponential scale
%     prob = tabulate(r);
%     nonzeop = prob(prob(:,3) > 0,3)/100;
%     h = sum(nonzeop.*log2(nonzeop));

%% in log scale
%     nr = 3;    %number of receptor
%     nl = 1;  %number of odor
%     EPS = 1e-3;
%     base = 4; % discret response
%     wnew = reshape(w,nr,nl);
%     response = zeros(nr,max(size(trainData)));
%     h = 0;
%     for i0 = 1:nr
%         exci = exp(wnew(i0,:))*trainData;
%         response(i0,:) = exci./(1+exci);
%         h = h - HShannon_KDP_estimation(response(i0,:),col);
%     end
%     %h = - HShannon_KDP_estimation(response,col);
% %     penality = -0.5*sum(abs(w - mean(w)));
% %     h = h + penality;
% 
%%  define effective dimension as regularization
if param.regul
    alpha = 1.0/param.nRecep; %penality of single channel
    C = cov(resp');
    h = h -alpha*trace(C)^2/trace(C*C);
end

if param.regularize
    lambda = 1;
    C = lambda*max(((max(max(w)))-5*param.lSig),0)*sum(sum(w(w>0).^2));
    h = h + C;
%% regularization, distance from uniform distribution
% if regul
%     distUnif = 0;
%     for i0 = 1:size(resp,1)
%         [Y,T] = ecdf(resp(i0,:));
%         distUnif = distUnif + max(abs(Y-T))/size(resp,1);
%     end
%     h  = h + distUnif;
% end

end

function h = entropyStructure(w)
% test the effect of maximizing entropy on structured data
global trainData col param
resp = zeros(2,1e4);
w = reshape(w,2,3);
if strcmp(param.function,'mm')
   resp = 2*w*trainData./(w*ones(size(trainData))+w*trainData);
   resp(resp==1) = 1-abs(normrnd(0,1e-10,size(find(resp==1))));
else
    resp = w*trainData./(max(max(w*trainData)));
     resp(resp==1) = 1-abs(normrnd(0,1e-10,size(find(resp==1))));
      resp(resp==0) = abs(normrnd(0,1e-10,size(find(resp==0))));
     
end
h = -HShannon_KDP_estimation(resp,col);
end


function h = entropyCostExciInhi(x)
% this function return differential entropy 
% interaction matrix is sampled from x, typcially the mean and std
% both excitation and inhibition has been consiered
% where the vecoterized weight matrix elements are reshaped into two matrix
% corresponding to excitation and hibition
% the interaction function contain parameter alpha, which modulate the base
% line activity

% define global variables
global trainData  param
persistent  col mult
mult = 1;
col = HShannon_KDP_initialization(mult);

% reshap x into two matrix
we = reshape(x(1:param.nOdor*param.nRecep),[param.nRecep,param.nOdor]);   %excitation 
wi = reshape(x((1+param.nOdor*param.nRecep):end),[param.nRecep,param.nOdor]); % inhibition

% response to training stimuli
resp = respExciInhiFun(trainData,we,wi,param.alpha);
h =  - HShannon_KDP_estimation(resp,col);

% regularization
% f = -h/(1+h)*sum(abs(we(:) - wi(:)))/param.nOdor/param.nRecep;
f = -sum(abs(we(:) - wi(:)))/param.nOdor/param.nRecep;
h = f + h;

end
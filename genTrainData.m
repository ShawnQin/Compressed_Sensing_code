function trainData = genTrainData(param,varargin)
% this function generate lognormal distributed odor mixture
% param the data structure with all the parameter shared by main program
% varargin  extral parameter to specify the correlation matrix

%%
if param.sparsity
    spars = ceil(param.spRatio*param.nOdor); %number of none zero elements
else
    spars = param.nOdor;
end

if param.corr
   if nargin == 1
      if isempty(param.eig)
          eigVal = specifyEig(param.nOdor); %random eigen values
      else
          eigVal = specifyEig(param.nOdor,param.eig);
      end
      corrCoefMat = randCorrCoef('buildin',eigVal);
   elseif nargin ==2
       corrCoefMat = varargin{1};
   elseif nargin > 2
       corrType = varargin{2};
       pm = varargin{2};
       corrCoefMat = randCorrCoef(corrType,pm);
       
   end
   
else
   corrCoefMat = diag(ones(1,param.nOdor)); %without correlation
%            trainData = abs(normrnd(mu,sig,[numOdor,numSamp]));
end

if strcmp(param.dType, 'binary')
   nSamp = min(param.nSamp,nchoosek(param.nOdor,spars));
   trainData = generateTrainData(nSamp,param.nOdor,corrCoefMat,...
               [],[],param.dType,spars,param.spType);
else
    if strcmp(param.dType, 'gauss')
        sig = param.gSig;
        mu = param.gMu;
    elseif strcmp(param.dType, 'lognorm')
        sig = param.lSig;
        mu = param.lMu;
    elseif strcmp(param.dType, 'powerlaw')
        sig = param.lSig; % the exponent of power law
        mu  = param.c0;  % the minimum concentration 
    elseif strcmp(param.dType, 'exponential')
        sig = param.lSig; % sig is just the lambda in exponential distribution
        mu = 0;  %in this case, only one parameter is needed
    end
    trainData = generateTrainData(param.nSamp,param.nOdor,corrCoefMat,...
               mu,sig,param.dType,spars,param.spType);
end
    


% if strcmp(param.dType,'gauss')
%        if param.corr
%            if isempty(param.eig)
%                eigVal = specifyEig(param.nOdor,param.nRecep); %random eigen values
%            else
%                eigVal = specifyEig(param.nOdor,param.nRecep,param.eig);
%            end
%            corrCoefMat = randCorrCoef('buildin',eigVal);
% %            trainData = generateTrainData(param.nSamp,param.nOdor,corrCoefMat,...
% %                param.gMu,param.gSig,'gauss',spars);
%        else
%            corrCoefMat = diag(ones(1,param.nOdor)); %without correlation
% %            trainData = abs(normrnd(mu,sig,[numOdor,numSamp]));
%        end
%        trainData = generateTrainData(param.nSamp,param.nOdor,corrCoefMat,...
%                param.gMu,param.gSig,'gauss',spars);
% elseif strcmp(dataType,'lognorm')
%         mu = 0;
%         sig = 2; 
%         % with correlation, although the method should be examined
%         if correlation
%             %eigVal = specifyEig(numOdor,numRecep,[7,5,4,3]);  % fixed the first 4 dimension for repeatable
%             eigVal = specifyEig(numOdor,numRecep);  % fixed the first 4 dimension for repeatable
%             corrCoefMat = randCorrCoef('buildin',eigVal);
%             trainData = generateTrainData(numSamp,numOdor,corrCoefMat,mu,sig,'lognorm'); %test without sparse
%         % default, without correlation
%         else
%             inx = datasample(1:numOdor*numSamp,numSamp*spars,'Replace',false);
%             trainData(inx) = exp(normrnd(mu,sig,[numSamp*spars,1]));
%         end
% end
% %%
% %trainData = zeros(nDim,numSamp);
% SIG = param.sig*ones(1,param.nOdor);
% MU = mu*ones(1,nDim);
% VAR = corrCoefMat.*(SIG'*SIG);
% if isempty(varargin)
%     trainData = exp(mvnrnd(MU,VAR,numSamp));    
% else
%     %sparsity
%     sp = varargin{1};  %number of nonzero elements
%     if sp >= nDim
%         error('sparsity mast be smaller than nDim!')
%     else
%         trainData = zeros(nDim,numSamp);
%         nonZeroInx = datasample(1:numSamp*nDim,numSamp*sp); %non zero index
%         temp = exp(mvnrnd(MU,VAR,numSamp));
%         trainData(nonZeroInx) = temp(nonZeroInx);     
%     end
% end
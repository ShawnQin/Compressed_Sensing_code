function trainData = generateTrainData(numSamp,nDim,corrCoefMat,mu,sig,dType,sp,spType,varargin)
% this function generate lognormal distributed odor mixture
% numSamp  number of samples, default 1e4
% nDim     dimension of odor space
% corrCoefMat  correlation coefficient matrix 
% mu       average of logrithm concentration, default 0
% sig      std of logrithm concentration, default 2
% varargin  extral parameter to specify sparsity

%trainData = zeros(nDim,numSamp);
% SIG = sig*ones(1,nDim);
% MU = mu*ones(1,nDim);
% VAR = corrCoefMat.*(SIG'*SIG);

% note: for power law distribution, sig is just the exponent
%       and mu is just the minium concentration


if length(mu) > 1
    if length(mu) == nDim && length(sig) == nDim
        MU = reshape(mu,1,nDim);
        SIG = reshape(sig,1,nDim);
    else
        error('dimension of mu, sig donot match with nDim!')
    end
else
    if length(sig) == 1
        SIG = sig*ones(1,nDim);
        MU = mu*ones(1,nDim);
    else
        %this is actually unecessary
        SIG = 1;
        MU = 0;
    end
end
VAR = corrCoefMat.*(SIG'*SIG);

if sp == nDim %sparse input
    if strcmp(dType,'lognorm')
        trainData = transpose(exp(mvnrnd(MU,VAR,numSamp)));
    elseif strcmp(dType,'gauss')
        trainData = transpose(abs(mvnrnd(MU,VAR,numSamp)));
    elseif strcmp(dType,'powerlaw')
        temp = rand(nDim,numSamp);
        trainData = mu*(1-temp).^(1/(1-sig)); % here sig is just the exponent of power law
    elseif strcmp(dType,'exponential')
        trainData = exprnd(sig,nDim,numSamp);
    end
        
else
    %sparsity
%     sp = varargin{1};  %number of nonzero elements
    if sp >= nDim
        error('sparsity mast be smaller than nDim!')
    else
        trainData = zeros(nDim,numSamp);
%         nonZeroInx = datasample(1:numSamp*nDim,numSamp*sp); %non zero index
        if strcmp(dType,'lognorm')
            temp = transpose(exp(mvnrnd(MU,VAR,numSamp)));
        elseif strcmp(dType,'gauss')
            temp = transpose(abs(mvnrnd(MU,VAR,numSamp)));
        elseif strcmp(dType,'powerlaw')
            temp = mu*(1 - rand(nDim,numSamp)).^(1/(1-sig));
        elseif strcmp(dType,'exponential')
            temp = exprnd(sig,nDim,numSamp);
        elseif strcmp(dType,'binary')
            
            if numSamp < nchoosek(nDim,sp)
                temp = ones(nDim,numSamp);
            else
                numSamp = nchoosek(nDim,sp);
                inxM = nchoosek(1:nDim,sp);
                trainData = zeros(nDim,numSamp);
                for i0 = 1:numSamp
%                     ix = inxM(i0,:);
                    trainData(inxM(i0,:),i0) = 1;
                end
                return
            end
        end
        
        if strcmp(spType,'absolute')  %absolute number
            for i0 = 1:numSamp
                indx = randperm(nDim,sp);
                trainData(indx,i0) = temp(indx,i0); 
            end 
        elseif strcmp(spType,'average')  %average number
            indx = randperm(numSamp*nDim,numSamp*(nDim-sp));
            temp(indx) = 0;
            trainData = temp;
        end

    end
end
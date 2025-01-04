function trainData = trainDataGen(numSamp,N,sp,sig)
%     MU = zeros(N,1);
%     SIG = sig*ones(1,N);
%     VAR = SIG'*SIG;
    temp = exp(sig*randn(N,numSamp));
    trainData = zeros(N,numSamp);
    for i0 = 1:numSamp
        indx = randperm(N,sp);
        trainData(indx,i0) = temp(indx,i0); 
    end 
%     indx = randperm(numSamp*N,numSamp*(N-sp));
%     temp(indx) = 0;
%     trainData = temp;
end
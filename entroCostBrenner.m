function h = entroCostBrenner(s)
    
% define persistent varible to save time
    global trainData param
    nr = param.nRecep;
    nl = param.nOdor;

    if strcmp(param.method,'wholeMat')
    
        if param.noise   %adding noise, default Gaussian
           r = int8(reshape(s,[nr,nl])*trainData + normrnd(0,param.noiseSig,nr,param.nSamp) >= 1); 
        else
           r = int8(reshape(s,[nr,nl])*trainData >= 1);
        end
              
        pt = tabulate(bi2de(r')); % a vector of different pattern in decimal number
        nonzeop = pt(pt(:,3) > 0,3)/100;
        h = sum(nonzeop.*log2(nonzeop));
        
        % whether with regularization or not
        if param.regularize
            lambda = 1;
            C = lambda*max((max(s)-5*param.lSig),0)*sum(s(s>0).^2);
            h = h + C;
        end
        
        
    elseif strcmp(param.method,'dist')
        N = 20; %entropy is esimated as the average of ten random matrix
        if strcmp(param.dType,'binary')  
            %sparse binary matrix
            h = 0;
            for i0 = 1:N
                W = zeros(nr,nl);
                W(randperm(nr*nl,round(nr*nl*s))) = 1;
                r = int8(W*trainData>= 1);
                pt = tabulate(bi2de(r')); % a vector of different pattern in decimal number
                nonzeop = pt(pt(:,3) > 0,3)/100;
                h = h + sum(nonzeop.*log2(nonzeop))/N; 
            end
        elseif strcmp(param.dType,'lognorm')
            % lognormal distributed
            h = 0;
            for i0 = 1:N
                W = zeros(nr,nl);
                W(randperm(nr*nl,round(nr*nl*s(3)))) = exp(normrnd(s(1),s(2),round(nr*nl*s(3)),1));
                r = int8(W*trainData >= 1);
                pt = tabulate(bi2de(r')); % a vector of different pattern in decimal number
                nonzeop = pt(pt(:,3) > 0,3)/100;
                h = h + sum(nonzeop.*log2(nonzeop))/N; 
            end
               
        end
    end
        
end
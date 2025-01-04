function S = randCorrCoef(varargin)
% random correlation coefficient generator
% combine several methods together 
% vine, onion and factor method
if isempty(varargin)
    error('You have to specify the method: vine,onion or factor!')
else
    if ~ischar(varargin{1})
        error('the first argument has to be a char!')
    elseif(strcmp(varargin{1},'vine'))
        % using vine method,
        % d  is the dimension
    % eta   the parameter in beta distribution, equal alpha and beta
    %       larger eta correspond to stronger correlation
        if nargin < 3
            error('two real numbers are needed for this method work!')
        else
            d = varargin{2};
            eta = varargin{3};
            beta = eta + (d-1)/2;   
            P = zeros(d);           %// storing partial correlations
            S = eye(d);

        for k = 1:d-1
            beta = beta - 1/2;
            for i = k+1:d
                P(k,i) = betarnd(beta,beta); %// sampling from beta
                P(k,i) = (P(k,i)-0.5)*2;     %// linearly shifting to [-1, 1]
                p = P(k,i);
                for l = (k-1):-1:1 %// converting partial correlation to raw correlation
                    p = p * sqrt((1-P(l,i)^2)*(1-P(l,k)^2)) + P(l,i)*P(l,k);
                end
                S(k,i) = p;
                S(i,k) = p;
            end
        end
        %// permuting the variables to make the distribution permutation-invariant
        permutation = randperm(d);
        S = S(permutation, permutation);
        end
    elseif(strcmp(varargin{1},'onion'))
        %// ONION METHOD to generate random correlation matrices distributed randomly
        d = varargin{2};
        S = 1;       
        for k = 2:d
            y = betarnd((k-1)/2, (d-k)/2); %// sampling from beta distribution
            r = sqrt(y);
            theta = randn(k-1,1);
            theta = theta/norm(theta);
            w = r*theta;
            [U,E] = eig(S);
            R = U*E.^(1/2)*U';             %// R is a square root of S
            q = R*w;
            S = [S q; q' 1];               %// increasing the matrix size
        end
    elseif(strcmp(varargin{1},'factor'))
        % d is the dimension of covairance coefficient matrix
        % k is interger equal or smaller than d, which control the elements of the
        %   matrix, larger k correspond to weak correlation
        if nargin < 3
            error('two real numbers are needed for this method work!')
        else
            d = varargin{2};
            k = varargin{3};
            W = randn(d,k);
            S = W*W' + diag(rand(1,d));
            S = diag(1./sqrt(diag(S))) * S * diag(1./sqrt(diag(S)));
        end
        
    elseif(strcmp(varargin{1},'buildin'))
        %using build in function 'gallery' to generate correlation matrix
        eigValue = varargin{2};  %dimension of correlation matrix
        S = gallery('randcorr',eigValue);
    elseif(strcmp(varargin{1},'exp'))
        % specify the eigenvalues as exponentially distributed, this method
        % is proposed in the paper Tkacik et al, PNAS,2010
        
        d = varargin{2}; %dimension of the matrix
        if nargin > 2
            lbd = varargin{3};  %specify how the change of exponential eigenvalues
            if lbd <=0
                error('lambda has to be a positive number!')
            end
        else
            % default lambda of exponential distribution
            lbd = 1/log(2);
        end
        %first generate exponentially distributed eigenvalue of
        %correlationmatrix       
        eigVal = exprnd(lbd,d,1);
        
        %generate a matrix to get the transformation matrix
        temp = randn(d,d);
        A = (temp + temp')/2;
        [P,~] = eig(A);
        
        %get the covariance matrix
        C = P*diag(eigVal)*P';
        S = zeros(d,d);
        for i0 = 1:d
            for j0 = 1:d
                S(i0,j0) = C(i0,j0)/sqrt(C(i0,i0)*C(j0,j0));
            end
        end
        
    end
end
        
end
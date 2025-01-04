
function trainData = genSkewGauss(numSamp,nDim,sp,xi,omi,alpha,varargin)
% this problem generate skewed odor concentration distribution
% M,N    this dimension of the size
% sp,    sparsity of odor, an integer, should be larger than N
% xi,    the position of skewed Gaussian
% omi,   the shape parameter
% alp,   determine the skewness

% check the spareless parameter
if sp > nDim
    error('exceed the total number of possible odorants!')
end

% define the pdf of skewed gaussian distribution
gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
skewStand = @(x,alpha) 2*gaussian(x).*normcdf(alpha*x);

% skewedgaussian = @(x,xi,w,alpha) 2/w*gaussian((x-xi)/w).*normcdf(alpha*(x-xi)/w);
% generate total nonzero elements
temp = zeros(sp*numSamp,1);

count = 1;
Const = 2;  %rejection method
while count <= sp*numSamp
    u = rand;
    v = randn;
    if Const*u*gaussian(v) <=  skewStand(v,alpha)
        temp(count) = exp(omi*(v - alpha*sqrt(2)/sqrt(pi*(1+alpha^2))));  %shift to 0
        count = count + 1;
    end
end
temp = reshape(temp,sp,numSamp);
% generate sparse odor samples as
trainData = zeros(nDim,numSamp);

for i0 = 1:numSamp
    indx = randperm(nDim,sp);
    trainData(indx,i0) = temp(:,i0);
end


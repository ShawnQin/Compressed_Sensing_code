function mse = LostFunFit(p,varargin)
% if x is a vecotor two elements, use the mse
% if x is a scalar use the histogram method
% x(1)   p
% x(2)   c0
global newData

regul = true;


if length(p) > 1
    if regul
        lamd = 0.3; %penelty
        mse = sum((-(p(2) - abs(newData(:,1)).^p(1)).^(1/p(1)) - newData(:,2)).^2) + lamd*abs(p(1));
%     mse = (-(p(2)^p(1) - abs(newData(:,1)).^p(1)).^(1/p(1)) - newData(:,2)).^2;
    else
        mse = sum((-(p(2) - abs(newData(:,1)).^p(1)).^(1/p(1)) - newData(:,2)).^2);
    end
    
else
    mse = std((abs(newData(:,1)).^p + abs(newData(:,2)).^p).^(1/p));

end
end
function [f,vfGrad] = targetFun(w,X,Y)
% f is the returned target function
% X is the reference concentration of odorants
% Y is the observed Or input
% note that X and Y are in linear scale
wm = reshape(w,[],size(Y,1)+1);
Xhat = wm(:,1:end-1)*Y + wm(:,end)*ones(1,length(Y));
f = mean(sum((log(1+exp(Xhat)) - log(1+X)).^2));

% return the gradient
N = size(X,1);
S = size(X,2);
Yhat = 1./(1+exp(-Xhat));
temp = [(log(1+exp(Xhat)) - log(1+X)).*Yhat*Y',sum(log(1+exp(Xhat)) - log(1+X).*Yhat,2)]/N/S;
vfGrad = temp(:);
end
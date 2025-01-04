function entr = entropycostAnav2(w) 
global param
%using numerical integral for cost function, for 2 odor case.
% param.lSig =2;
% param.h =2;
% sig = 2;
% h =  1;
% J = @(y) reshape(reshape(exp(h*(w+reshape(kron(y,ones(2,1)),4,size(y,2)))),2,2*size(y,2))./(1+kron(((exp(reshape(w,2,2))*exp(y))).^h,ones(1,2))).^2,4,size(y,2));
% P = @(x) prod(normpdf(x,0,param.lSig));
% C = @(J) J(1,:).*J(4,:)-J(2,:).*J(3,:);
% fun = @(c1,c2) P([c1;c2]).*log(abs(C(J([c1;c2]))));
if param.nOdor == 1
   sig = param.lSig;  
   intFun = @(x) normpdf(x,0,sig).*log(abs(Jac1(w,x,param.h)));
   entr = -1/2*integral(intFun,-5*sig,5*sig);
elseif param.nOdor ==2
   sig1 = param.lSig;
   sig2 = param.lSig;
   intFun = @(x,y) normpdf(x,0,sig1).*normpdf(y,0,sig2).*log(abs(JacNew(w,x,y,param.h))+1e-10);
   entr = -1/2*integral2(intFun,-5*sig1,5*sig1,-5*sig2,5*sig2);
else
    error('odor number has too be 1 or two!')
end

% add regularization
if param.regularize
   lambda = 1;
   C = lambda*max((max(w(:))-10*param.lSig),0)*sum(sum(w(w>0).^2));
   entr = entr + C;
end

% intFun = @(x,y) normpdf(x,0,sig1).*normpdf(y,0,sig2).*log(abs(Jac1(w,x,y,param.h)));
% intFun = @(x) normpdf(x,0,sig1).*log(abs(Jac1(w,x,param.h)));

% tic
% for i0 = 1:10
%     h1 = -1/2*integral2(intFun,-5*sig,5*sig,-5*sig,5*sig);
% end
% toc
% entr = -1/2*integral(intFun,-5*sig1,5*sig1);
% entr = -1/2*integral2(intFun,-5*sig1,5*sig1,-5*sig2,5*sig2);
%  f = chebfun2(fun,-3,3,-3,3,'vectorize');
%  h = -sum2(f);
% h = -quad2d(fun,-10,10,-10,10);
%   a = 0.1;
%   c = -5*param.lSig:a:5*param.lSig;%set the range of grid points.
%  [c1,c2] = ndgrid(c);
%   c = [c1(:),c2(:)]';
% %  h = 0;
% % for i = 1:length(c)
% tic
% for i0  = 1:10
%  h = -a^2*sum(fun(c(1,:),c(2,:)));
%  disp(h)
% end
% toc
%  h = h+a^2*real(h0);
% end%numerically calculating the integral.
% h = -h;
% disp(entr)
%% direct integration using matlab integral2

end

function J = Jac(W,x,y,h)
w = reshape(W,2,2);
% J = zeros(2,2);
t1 = exp(h*(w(1,1) + x))./(1+exp(h*(w(1,1)+x))+exp(h*(w(1,2)+y))).^2;
t2 = exp(h*(w(1,2) + y))./(1+exp(h*(w(1,1)+x))+exp(h*(w(1,2)+y))).^2;
t3 = exp(h*(w(2,1) + x))./(1+exp(h*(w(2,1)+x))+exp(h*(w(2,2)+y))).^2;
t4 = exp(h*(w(2,2) + y))./(1+exp(h*(w(2,1)+x))+exp(h*(w(2,2)+y))).^2;
J = t1.*t4 - t3.*t2;
end

% when N = 2, M>=2
function chi = JacNew(W,x,y,h)
[row,col] = size(x);  %grid of X
nr = length(W)/2;     %number of receptors
xvec = x(:);          %vectorize grid of x
yvec = y(:);          %vectorize grid of y
newVal = [xvec';yvec'];

w = reshape(W,[],2);  %reshape W into a matrix


% nn = exp(h*(W+kron([xvec';yvec'],ones(nr,1))));
% dn = (1+reshape(kron(((exp(w)*exp(newVal))).^h,ones(1,2)),2*nr,[])).^2;
% J = nn./dn;    %Jacobian
J = exp(h*(W+kron([xvec';yvec'],ones(nr,1))))./(1+reshape(kron(((exp(w)*exp(newVal))).^h,ones(1,2)),2*nr,[])).^2;


% determinant
temp = sum(J(1:nr,:).^2,1).*sum(J((nr+1):end,:).^2,1) - (sum((J(1:nr,:).*J((nr+1):end,:)),1)).^2;
chi = reshape(temp,row,col);   %reshape back to the size of x or y

end

% when N=1, M>=2
function chi=Jac1(w,x,h)
    nr = length(w);
    temp = exp(h*(w+kron(x,ones(nr,1))))./(1+exp(h*(w+kron(x,ones(nr,1))))).^2;
    chi = sum(temp.^2,1);
%     chi = sum(temp.^2,1) + (0.01/2)^2;  % debug on 02/28/2018
end
function val = respFund1(W,X)
% this function return the Jacobian matrix of repsonse function
% W     an MxN matrix, log scale
% X     an N dimensional vector,linear scale
% val   an NxN matrix

% last modified 02/28/2018
% ==================================================================

% check dimension consistance
[row, col] = size(W);
if length(X) ~= col
    error('input dimensions of W and X are not consistent!')  
end

%  the following calculation uses the broadcast function of matrix
%  manipulation in MATLAB (automatic dimension match)
denomenator = (1+exp(W)*X).^2;
J = exp(W).*X'./denomenator;
val = J'*J;
  
end
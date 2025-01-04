function v = digitFun(X,n,h)
% X  a vector or matrix
% n  how many digit you want to convert
% h  control the slop of sigmoid
T = 1./(1+X.^h);
v = floor(T*n);
end
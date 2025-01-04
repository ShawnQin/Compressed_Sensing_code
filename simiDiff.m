function df = simiDiff(A,S0,a,varargin)
% this function return the difference between W and reference simiarity
% matirx

% A       M x N matrix
% S0       N x N matrix
% a       penalty weight 
% df      a weight difference between actual similarity matrix and
%          reference  matrix
N = size(A,2);  % number of odors

% calculate the actual "similarity matrix" 
T = sqrt(sum(A.*A,1)); %length of each vector
S = zeros(N,N);
df = 0;
for i0 =  1:(N-1)
    for j0  = (i0 + 1):N
        S(i0,j0) = dot(A(:,i0),A(:,j0))/T(i0)/T(j0);
        temp = S(i0,j0) - S0(i0,j0);
        if temp < 0
            df = df - a*temp;
        end
    end
end
end
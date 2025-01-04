function eigVal = specifyEig(n,peig,varargin)
% this function specify the eigenvalues of a correlation matrix
% depends on "specifyEig.m"
% n the dimension of correlation matrix
% peig  principle eigenvalues, can be a integer or a vector

if nargin == 1
    temp = rand(n,1);
    eigVal = n*temp/sum(temp);
else    
    if length(peig) == 1 && isinteger(peig)
        m = peig; %first m eigenvalues
        temp1 = rand(m,1);
        eigVal(1:m) = 0.95*n*temp1/sum(temp1);
        temp2 = rand(n-m,1);
        eigVal(m+1:n) = 0.05*n*temp2/sum(temp2);        
    elseif length(peig)>1 || ~isinteger(peig)
        m = length(peig);
        if sum(peig)>=n
            error('sum of eigenvalues must be smaller than n!')
        else
            eigVal(1:m) = peig;
            remain = n - sum(peig);
            temp3 = rand(n-m,1);
            eigVal(m+1:n) = remain*temp3/sum(temp3);
        end

    end
    


    %the first m eigvalues are more than 95% of the total value
%     eigVal = zeros(n,1);
%     if ~isempty(varargin)
%         m = length(varargin{1});
%         %sepcify the first m eigenvalues
%         eigVal(1:m) = varargin{1};
%         remain = n - sum(eigVal(1:m));
%         temp3 = remain*rand(n-m,1);
% %         temp3 = sort([0;remain*rand(n-m-1,1);remain]);
%     else   
%         temp1 = 0.95*n*rand(m-1,1);
%         temp2 = sort([0;temp1;0.95*n]);
%         eigVal(1:m) = diff(temp2);
%         temp3 = sort([0;0.05*n*rand(n-m-1,1);0.05*n]);
%     end
%     
%     
% %     eigVal(m+1:n) = diff(temp3);
%     eigVal(m+1:n) = temp3/sum(temp3);
    eigVal = sort(eigVal,'descend');
end

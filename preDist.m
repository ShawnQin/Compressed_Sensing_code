function D = preDist(N,varargin)
% this function returns a distance matrix
% there are two metric of distance, Euclidian or cosine

% N             number of odorants
% varargin{1}   distance metric, euclidian, cosine, and pearson
% varargin{2}   type of distance distribution, default beta
% varargin{3}   select partial or whole, default: whole

% D is a N by N symetric matrix
% distance metric type
if nargin > 1
    metric = varargin{1};
% else
%     metric = 'euclidian';
end

% draw the distance from a certian distribution
if nargin > 2
    distrType = varargin{2};
else
    distrType = 'beta';
end

if nargin > 4
    a = varargin{4};
    b = varargin{5};
else
    if strcmp(distrType,'beta')
        a = 3;
        b = 3;
    end
end

% select part all the whole matrix
if nargin > 3
    NS = varargin{3};
    if ~isnumeric(NS)
        error('selected index has to be integers, a vector')
    end
else
    NS = N;  % select the whole matrix
end

% Euclidian distance
if strcmp(metric,'euclidian')
    ranges = [0,1]; % the distance range    
elseif strcmp(metric,'cosine') || strcmp(metric,'pearson')
    ranges = [-1,1];   
else
    error('Unsupported distance metric! Metric has to be euclidian, cosine or pearson!')
end

% distribution of pairwise distance
% considering the first slected
D = zeros(NS,NS);
if strcmp(distrType,'beta')
%     Y = (ranges(2) - ranges(1))*(betarnd(a,b,N,N) - 0.5); %scale to [-1,1]
    for i0 = 1:(NS-1)
        for j0 = i0:NS
%             D(i0,j0) =  (ranges(2) - ranges(1))*(betarnd(a,b) - 0.5);
            D(i0,j0) = (-1)^binornd(1,0.5)*(betarnd(a,b)/2 + 0.5);  %revised on 6/4/2018
            D(j0,i0) = D(i0,j0);
        end
    end
    D(logical(eye(NS))) = 1;
elseif strcmp(distrType,'gauss')
    % trucated
elseif strcmp(distrType,'rand')
    for i0 = 1:(NS-1)
        for j0 = i0:NS
            D(i0,j0) =  (ranges(2) - ranges(1))*(rand - 0.5);
            D(j0,i0) = D(i0,j0);
        end
    end
    D(logical(eye(NS))) = 1;
elseif strcmp(distrType,'block')
    % the correlation inside blocks
    D = zeros(N,N);
    NB = max(4,floor(N/5));  % number of blocks
    BS = floor(N/NB);  % size of block
    
    %also use beta distribution
    a = 2;
    b = 1;
    for i0 = 1:NB
        for j0 = (BS*(i0 - 1) + 1):(BS*i0-1)
            for k0 = j0+1:BS*i0
                D(j0,k0) = betarnd(a,b);
            end
        end
    end
 elseif strcmp(distrType,'gradient')  %added on 06/08/2018
     % correlation decay with index distance, this is used when comparing
     % with Guangwei's work
     D = zeros(N,N);
     NB = 2;  % two blocks
     BS = floor(N/NB);  % size of block
     for i0 = 1:NB
         for j0 = (BS*(i0 - 1) + 1):(BS*i0-1)
             for k0 = j0+1:BS*i0
                 D(j0,k0) = exp(-abs(j0 - k0)/5);
             end
         end
     end
         
     
end


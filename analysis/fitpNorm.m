function [p,c0] = fitpNorm(data)
% this function fit the correlated function of w1 and w2 using CMA-ES
% method. It first transform the data into w1^p + w2^p, and minimize the
% variation of these transformed data. You can also use least squre error
% as targe function.
% this a childran function of 'analyzeOptMat_v2.m'

% data     the raw data which has to be preprocessed inorder to select the
%          desired data point along -3\sigma ~ 0

% last revised on 4/13/2018



%% first determine how many ring structure this data have
global newData

% this can be done manually (interactively)

% scatter plot
if size(data,2) ~= 2
    error('the dimesion of data do not matach!')
else
%     thrd = -3;   % in unit of sigma_c
%     temp = data(all(data > thrd,2),:);   %normalized to sigma_c
%     scatter(temp(:,1),temp(:,2));
    scatter(data(:,1),data(:,2));
    xlim([-5,1])
    ylim([-5,1])
    prompt = 'Input the number of band structure, an integer:';
    n = input(prompt);
    
    % initial guess of p
    p0 = 0.6;
    X = -4:0.02:0;
    %cycle through different band
    sepSeq = zeros(n,1); %distance threshold used to sparate different bands
    p = zeros(n,1);
    c0 = zeros(n,1);
    for i0 = 1:n
%         thrd = -2 - i0;  %threshold depends on the number of bands
        thrd = -3;       % modified on 4/18/2018
        %check if the band can be separated
        flag = 0;
        while flag == 0
            ix1 = all(data > thrd,2);
            ix2 = all(data < 3,2);  %get rid of too large value
            temp = data(all([ix1,ix2],2),:);  % last revised on 4/18
%             temp = data(all(data > thrd,2),:);
            figure(1)
            scatter(temp(:,1),temp(:,2));
            hold on
            prompt = 'Input the distance to separte the first one with others:';
            ds = input(prompt);
            Y = -(ds - abs(X).^p0).^(1/p0);
            plot(X,Y,'r-')
            prompt = 'Are the separatable now? y or n';
            fg = input(prompt,'s');
            if strcmp(fg,'y')
%                 hold off
                flag = 1;   %sucessed flag                
                sepSeq(i0) = ds;
                disp(n)
                if flag  == 1
                    if i0==1
                        inxSelect = abs(temp(:,1)).^p0  + abs(temp(:,2)).^p0 < sepSeq(i0);
%                     elseif i0==2
%                         inxSelect = temp(:,1).^p0  + temp(:,2).^p0 > sepSeq(1);
                    else
                        inx1 = all(abs(temp(:,1)).^p0  + abs(temp(:,2)).^p0 < sepSeq(i0),2);
                        inx2 = all(abs(temp(:,1)).^p0  + abs(temp(:,2)).^p0 > sepSeq(i0-1),2);
                        inxSelect = all([inx1,inx2],2);
                    end
                        
                    newData = temp(inxSelect,:);
                    
                    % using CMA-ES algorithm to find the best fitting
                    opts = cmaes;
%                     opts.LBounds = 0.01;
%                     opts.UBounds = 1;
                    opts.LBounds = [0;0];
                    opts.UBounds = [1;10];
                    opts.PopSize = 50;
                    x0 = [0.6,2.5];
%                     x0 = 0.6;
                    iniSig = [];
                    [xmin,fmin] = cmaes('LostFunFit',x0,iniSig,opts);
                    p(i0) = xmin(1);
                    c0(i0) = xmin(2);
%                     c0(i0) = mean(abs(newData(:,1)).^p(i0)  + abs(newData(:,2)).^p(i0));
                    
                    % plot the fitted curve on the scatter plot
                    X = -3:0.01:0;
                    Y = -(c0(i0) - abs(X).^p(i0)).^(1/p(i0));
                    plot(X,Y,'k-','LineWidth',3)
                    xlim([-3,0])
                    ylim([-3,0])
                    hold off
                    figure(1)
                    pause(1)   % just to allow you to check if fitting is good
                end
            
            
            end
             
        end
    end
    
end

end
% this program extract and merge all the data for each simulation run on
% the  linux cluster
%
close all
clear

% folder contain all th subfolders, default is the current folder
% folder = '/home/shan/Documents/MATLAB/theoNeurosci/olfaction/data/intN2_0411';
% folder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_Ndp';
folder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_skew';

allFNames = dir(folder);
folderFlag = [allFNames.isdir];
subFolders = allFNames(folderFlag);
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
dirNames = {subFolders.name};


% set the landmark in the folder names, so that we can just handle what we
% want. Right now we jut need the data folders that are from cluster
% N = 1;    %sparsity
% sig =2;   %standard deviation of

% using regular expression to find the relavent folders
str_marker = '2019-02-03';   %folder with this string contains the data we need
% str_marker = ['_S',num2str(S),'_sig',num2str(sig),'_2018-01-13']; 
% str_marker = ['_sig',num2str(sig),'_']; 

FIND = @(str) cellfun(@(c) ~isempty(c), regexp(dirNames,str,'once'));


%for the amplitude data
targetNames = dirNames(FIND(str_marker));

%extract and merge the mat file in the folders
%save file name
dateInfo = datestr(now,'yyyy-mm-dd');
for i0 = 1:length(targetNames)
    allFile = dir([folder,filesep,targetNames{i0},filesep,'*.mat']);
    files = {allFile.name}';
    
    if isempty(files)
        error('folder is empty!')
    else
        
        %cycle throught different N and R
        % get the valu of N and R
        s1 = '(?<= *N)[\d.]+(?=_)';
        s2 = '(?<= *_R)[\d]+(?=_)';
        s3 = '(?<= *_S)[\d]+(?=_)';
        s4 = '^\w{3,4}(?=N)';     
        s5 = '(?<= *sig)[\d.]+(?=_)';
        s6 = '(?<= *alp)[-\d.]+(?=_)';
%         s5 = '(?<= *lbd)[\d.]+(?=_)';  % for exponential input
        
        numOdor = str2num(char(regexp(targetNames{i0},s1,'match')));
        numRecp = str2num(char(regexp(targetNames{i0},s2,'match')));
        sp = str2num(char(regexp(targetNames{i0},s3,'match')));
        tp = char(regexp(targetNames{i0},s4,'match'));
        sig = str2num(char(regexp(targetNames{i0},s5,'match')));
        alp = str2num(char(regexp(targetNames{i0},s6,'match')));
        allMat = zeros(numOdor*numRecp,length(files));
        allfmin = zeros(length(files),1);
        
%         allCorrMat = cell(length(files),1); %store the correlation matrix
        
        for j0 = 1:length(files)
            temp = load(char(fullfile(folder,filesep,targetNames{i0},filesep,files{j0})));
            allMat(:,j0)  = temp.wmin;
            allfmin(j0) = temp.fmin;
%             allCorrMat{j0} = temp.corrCoefMat;
%             allCorrMat{j0} = temp.corrMat;
        end
    end
    

%       save([folder,filesep,'Ising_N',num2str(numOdor),'_R',num2str(numRecp),'_S',num2str(sp),'_sig',num2str(sig),...
%           '_',dateInfo,'.mat'],'allMat','allfmin','allCorrMat')
%       save([folder,filesep,'Gcmi_N',num2str(numOdor),'_R',num2str(numRecp),'_S',num2str(sp),'_sig',num2str(sig),...
%           '_',dateInfo,'.mat'],'allMat','allfmin')
      save([folder,filesep,'skew_N',num2str(numOdor),'_R',num2str(numRecp),'_S',num2str(sp),'_sig',num2str(sig),...
          '_alp',num2str(alp),'_',dateInfo,'.mat'],'allMat','allfmin')
end

%% extract when optimize the distribution
% this program extract and merge all the data for each simulation run on
% the  linux cluster
%
close all
clear

% folder contain all th subfolders, default is the current folder
folder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_distr';

allFNames = dir(folder);
folderFlag = [allFNames.isdir];
subFolders = allFNames(folderFlag);
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
dirNames = {subFolders.name};


% set the landmark in the folder names, so that we can just handle what we
% want. Right now we jut need the data folders that are from cluster

% using regular expression to find the relavent folders
str_marker = '2018-08-28';   %folder with this string contains the data we need

FIND = @(str) cellfun(@(c) ~isempty(c), regexp(dirNames,str,'once'));

%for the amplitude data
targetNames = dirNames(FIND(str_marker));

%extract and merge the mat file in the folders
%save file name
dateInfo = datestr(now,'yyyy-mm-dd');
for i0 = 1:length(targetNames)
    allFile = dir([folder,filesep,targetNames{i0},filesep,'*.mat']);
    files = {allFile.name}';
    
    if isempty(files)
        error('folder is empty!')
    else
        
        %cycle throught different N and R
        % get the valu of N and R
        s1 = '(?<= *N)[\d.]+(?=_)';
        s2 = '(?<= *_R)[\d]+(?=_)';
        s3 = '(?<= *_S)[\d]+(?=_)';
%         s4 = '^\w{3,4}(?=N)';
        s5 = '(?<= *sig)[\d.]+(?=_)';
%         s5 = '(?<= *lbd)[\d.]+(?=_)';  % for exponential input
        
        numOdor = str2num(char(regexp(targetNames{i0},s1,'match')));
        numRecp = str2num(char(regexp(targetNames{i0},s2,'match')));
        sp = str2num(char(regexp(targetNames{i0},s3,'match')));
%         tp = char(regexp(targetNames{i0},s4,'match'));
        sig = str2num(char(regexp(targetNames{i0},s5,'match')));
        allParam = zeros(3,length(files));
        allfmin = zeros(length(files),1);
        
        for j0 = 1:length(files)
            temp = load(char(fullfile(folder,filesep,targetNames{i0},filesep,files{j0})));
            allParam(:,j0)  = temp.wmin;
            allfmin(j0) = temp.fmin;
        end
    end
    
    save([folder,filesep,'gcmi_distr_N',num2str(numOdor),'_R',num2str(numRecp),'_S',num2str(sp),'_sig',num2str(sig),...
          '_',dateInfo,'.mat'],'allParam','allfmin')
end


%% extract and merge data for the one odor and many receptors

folder = './';
allFNames = dir(folder);
folderFlag = [allFNames.isdir];
subFolders = allFNames(folderFlag);
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
dirNames = {subFolders.name};

str_prefix = 'oneByM';   %folder with this string contains the data we need
% str_marker = ['_sig',num2str(sig),'_2018-01-14']; 
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(dirNames,str,'once'));


%for the amplitude data
targetNames = dirNames(FIND(str_prefix));
dateInfo = datestr(now,'yyyy-mm-dd');
for i0 = 1:length(targetNames)
    allFile = dir(fullfile(targetNames{i0},filesep,'*.mat'));
    files = {allFile.name}';
    if isempty(files)
        error('folder is empty!')
    else
        
        %cycle throught different N and R
        % get the valu of N and R
%         s1 = '(?<= *LogN)[\d.]+(?=_)';
%         s2 = '(?<= *_R)[\d]+(?=_)';
%         numOdor = str2num(char(regexp(targetNames{i0},s1,'match')));
%         numRecp = str2num(char(regexp(targetNames{i0},s2,'match')));
%         allMat = zeros(numOdor*numRecp,length(files));
%         allfmin = zeros(length(files),1);
        
        for j0 = 1:length(files)
            temp = load(char(fullfile(folder,filesep,targetNames{i0},filesep,files{j0})));
            allMat(:,j0)  = temp.wmin;
            allfmin(j0) = temp.fmin;
        end
    end
    
    
    
    save(['oneByM','_sig',num2str(sig),'_',dateInfo,'.mat'],'allMat','allfmin')
end

%% extract data for overcomplete situation (M>N)
% gradient-based algorithm
N = 5;
sig = 4;
repeats = 20;
% R = [N:20,30,40,50,100,200];
% R = [3,4,5,8,10,13,15,19,20,30,40,50,100,200];
R = [6:20,30,40,100,200];

folder = './gradInfoMax/gradN5_1e3_0204';
allFNames = dir(folder);
folderFlag = [allFNames.isdir];
subFolders = allFNames(folderFlag);
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
dirNames = {subFolders.name};

str_prefix = ['_sig',num2str(sig)];   %folder with this string contains the data we need
% str_marker = ['_sig',num2str(sig),'_2018-01-14']; 
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(dirNames,str,'once'));


%for the amplitude data
targetNames = dirNames(FIND(str_prefix));
dateInfo = datestr(now,'yyyy-mm-dd');
allMat = cell(length(R),1);
allFmin = cell(length(R),1);

% order of receptor
s1 = '(?<= *_R)[\d.]+(?=_)';
% orderInx = zeros(length(targetNames));
Recp = zeros(length(targetNames),1);
for i0 = 1:length(targetNames)
   Recp(i0) = str2num(char(regexp(targetNames{i0},s1,'match')));
end
[~,orderInx] = sort(Recp);

for i0 = 1:length(targetNames)
    dataFile = char(fullfile(folder,targetNames{orderInx(i0)},'infoMax.mat'));
    if exist(dataFile,'file')
        
    load(dataFile);
    
    % delete some "bad" data
    badInx = [];
    allFmin{i0} = [];
    for j0 = 1:repeats
        if length(totalCost{j0}) < 10 || isnan(totalCost{j0}(end))
            badInx = [badInx,j0]; %store bad data
        else
            allFmin{i0} = [allFmin{i0};totalCost{j0}(end)];
        end
%         if isinf(real(fmin(i0))) || isnan(fmin(i0))
%             badInx = [badInx,j0]; %store bad data
%         else
%             allFmin{i0} = [allFmin{i0};real(fmin(i0))];
%         end
       
    end
   allMat{i0} = allW; 
   if ~isempty(badInx)
        allMat{i0}(:,badInx) = [];
   end
   
    else
        allMat{i0} = nan(N*Recp(orderInx(i0)),repeats);
    end
    
end
save([num2str(N),'xM','_sig',num2str(sig),'_1e3_',dateInfo,'.mat'],'allMat','allFmin')

%% Extract data when both nhibition and excitation is present
% this program extract and merge all the data for each simulation run on
% the  linux cluster
%
close all
clear

% folder contain all th subfolders, default is the current folder
folder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_inhi/N50M10S2sig2_1013';
allFNames = dir(folder);
folderFlag = [allFNames.isdir];
subFolders = allFNames(folderFlag);
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
dirNames = {subFolders.name};


% set the landmark in the folder names, so that we can just handle what we
% want. Right now we jut need the data folders that are from cluster
% N = 1;    %sparsity
sig =2;   %standard deviation of

% using regular expression to find the relavent folders
str_marker = '2018-10-13';   %folder with this string contains the data we need
% str_marker = ['_S',num2str(S),'_sig',num2str(sig),'_2018-01-13']; 
% str_marker = ['_sig',num2str(sig),'_']; 

FIND = @(str) cellfun(@(c) ~isempty(c), regexp(dirNames,str,'once'));


%for the amplitude data
targetNames = dirNames(FIND(str_marker));

%extract and merge the mat file in the folders
%save file name
dateInfo = datestr(now,'yyyy-mm-dd');
for i0 = 1:length(targetNames)
    allFile = dir([folder,filesep,targetNames{i0},filesep,'*.mat']);
    files = {allFile.name}';
    
    if isempty(files)
        error('folder is empty!')
    else
        
        %cycle throught different N and R
        % get the valu of N and R
        s1 = '(?<= *N)[\d.]+(?=_)';
        s2 = '(?<= *_R)[\d]+(?=_)';
        s3 = '(?<= *_S)[\d]+(?=_)';
        s4 = '(?<= *_frac)[\d]+(?=_)';
%         s4 = '^\w{3,4}(?=N)';        
        s5 = '(?<= *sig)[\d.]+(?=_)';
        s6 = '(?<= *alp)[\d.]+(?=_)';
        numOdor = str2num(char(regexp(targetNames{i0},s1,'match')));
        numRecp = str2num(char(regexp(targetNames{i0},s2,'match')));
        sp = str2num(char(regexp(targetNames{i0},s3,'match')));
        frac = char(regexp(targetNames{i0},s4,'match'));
        sig = str2num(char(regexp(targetNames{i0},s5,'match')));
        alpha = str2num(char(regexp(targetNames{i0},s6,'match')));
        allMat = zeros(numOdor*numRecp,length(files));
        allSign = zeros(numOdor*numRecp,length(files));
        allfmin = zeros(length(files),1);
        
        for j0 = 1:length(files)
            temp = load(char(fullfile(folder,filesep,targetNames{i0},filesep,files{j0})));
            allMat(:,j0)  = temp.wmin;
            allSign(:,j0) = temp.Sign;

            allfmin(j0) = temp.fmin;
        end
    end
    
    
    
%     save(['u1m_N',num2str(numOdor),'_R',num2str(numRecp),...
%     '_S',num2str(S),'_sig',num2str(sig),'_',dateInfo,'.mat'],'allMat','allfmin')
%       save([folder,filesep,'spDepend_N',num2str(numOdor),'_R',num2str(numRecp),'_S',num2str(sp),'_sig',num2str(sig),...
%           '_',dateInfo,'.mat'],'allMat','allfmin')
      save([folder,filesep,'gcmiInhi_N',num2str(numOdor),'_R',num2str(numRecp),'_S',num2str(sp),'_sig',num2str(sig),...
          '_alp',num2str(alpha),'_frac',num2str(frac),'_',dateInfo,'.mat'],'allMat','allfmin','allSign')
end

%% Extract data with parital inhibition
% there are some fraction of receptors that can be both excited or
% inhibited by odorants
%
close all
clear

% folder contain all th subfolders, default is the current folder
folder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_partInhi';
allFNames = dir(folder);
folderFlag = [allFNames.isdir];
subFolders = allFNames(folderFlag);
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
dirNames = {subFolders.name};


% set the landmark in the folder names, so that we can just handle what we
% want. Right now we jut need the data folders that are from cluster

% using regular expression to find the relavent folders
str_marker = '2019-06-25';   %folder with this string contains the data we need

FIND = @(str) cellfun(@(c) ~isempty(c), regexp(dirNames,str,'once'));


%for the amplitude data
targetNames = dirNames(FIND(str_marker));

%extract and merge the mat file in the folders
%save file name
dateInfo = datestr(now,'yyyy-mm-dd');
for i0 = 1:length(targetNames)
    allFile = dir([folder,filesep,targetNames{i0},filesep,'*.mat']);
    files = {allFile.name}';
    
    if isempty(files)
        error('folder is empty!')
    else
        
        %cycle throught different N and R
        % get the valu of N and R
        s1 = '(?<= *N)[\d.]+(?=_)';
        s2 = '(?<= *_R)[\d]+(?=_)';
        s3 = '(?<= *_S)[\d]+(?=_)';
        s4 = '(?<= *frac)[\d.]+(?=_)';  % fraction of inhibiotry receptors
        s5 = '(?<= *sig)[\d.]+(?=_)';
        s6 = '(?<= *alp)[\d.]+(?=_)';
        numOdor = str2num(char(regexp(targetNames{i0},s1,'match')));
        numRecp = str2num(char(regexp(targetNames{i0},s2,'match')));
        sp = str2num(char(regexp(targetNames{i0},s3,'match')));
        frac = str2num(char(regexp(targetNames{i0},s4,'match')));
        sig = str2num(char(regexp(targetNames{i0},s5,'match')));
        alpha = str2num(char(regexp(targetNames{i0},s6,'match')));
%         allMat = zeros(numOdor*numRecp,length(files));
        allInhiW = zeros(numOdor*ceil(numRecp*frac),length(files));
        allExciW = zeros(numOdor*(numRecp - ceil(numRecp*frac)),length(files));
        allSign = zeros(numOdor*ceil(numRecp*frac),length(files));
        
        % store all the r0 if have, added 06/21/2019
        allr0 = zeros(ceil(numRecp*frac),length(files));
%         allInhiW = cell(length(files),1);
%         allExciW = cell(length(files),1);
%         allSign = cell(length(files),1);        
        allfmin = zeros(length(files),1);
        
        for j0 = 1:length(files)
            temp = load(char(fullfile(folder,filesep,targetNames{i0},filesep,files{j0})));
            allInhiW(:,j0)  = temp.wInhi(:);
            allExciW(:,j0) = temp.wp(:);
            allSign(:,j0) = temp.Sign(:);

            allfmin(j0) = temp.fmin;
            allr0(:,j0) = temp.allr0;  % added 06/21/2019
        end
    end
    
   
      save([folder,filesep,'gcmiPartInhi_N',num2str(numOdor),'_R',num2str(numRecp),'_S',num2str(sp),'_sig',num2str(sig),...
          '_alp',num2str(alpha),'_frac',num2str(frac),'_',dateInfo,'.mat'],'allInhiW','allExciW','allfmin','allSign')
end

%% extract and merge data of both excitation and inihibition from Qianyi's simulation
close all
clear

% folder contain all th subfolders, default is the current folder
%folder = '/Users/shan/Dropbox/olfactionProject/data/GcmiInhibit-06-13';
folder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_inhi';
allFNames = dir(folder);
folderFlag = [allFNames.isdir];
subFolders = allFNames(folderFlag);
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
dirNames = {subFolders.name};


% set the landmark in the folder names, so that we can just handle what we
% want. Right now we jut need the data folders that are from cluster
% N = 1;    %sparsity
% sig =2;   %standard deviation of

% using regular expression to find the relavent folders
str_marker = '2018-10-07';   %folder with this string contains the data we need
% str_marker = ['_S',num2str(S),'_sig',num2str(sig),'_2018-01-13']; 
% str_marker = ['_sig',num2str(sig),'_']; 

FIND = @(str) cellfun(@(c) ~isempty(c), regexp(dirNames,str,'once'));


%for the amplitude data
targetNames = dirNames(FIND(str_marker));

%extract and merge the mat file in the folders
%save file name
dateInfo = datestr(now,'yyyy-mm-dd');
for i0 = 1:length(targetNames)
    allFile = dir([folder,filesep,targetNames{i0},filesep,'*.mat']);
    files = {allFile.name}';
    
    if isempty(files)
        error('folder is empty!')
    else
        
        %cycle throught different N and R
        % get the valu of N and R
        s1 = '(?<= *N)[\d.]+(?=_)';
        s2 = '(?<= *_R)[\d]+(?=_)';
        s3 = '(?<= *_S)[\d]+(?=_)';
        s4 = '^\w{3,4}(?=N)';
        s5 = '(?<= *sig)[\d.]+(?=_)';
        s6 = '(?<= *alp)[\d.]+(?=_)';
        %s7 = '(?<= *R)[\d.]+(?=_2018)';
        numOdor = str2num(char(regexp(targetNames{i0},s1,'match')));
        numRecp = str2num(char(regexp(targetNames{i0},s2,'match')));
        sp = str2num(char(regexp(targetNames{i0},s3,'match')));
        tp = char(regexp(targetNames{i0},s4,'match'));
        sig = str2num(char(regexp(targetNames{i0},s5,'match')));
        alpha = str2num(char(regexp(targetNames{i0},s6,'match')));
        %r0 = str2num(char(regexp(targetNames{i0},s7,'match')));
        allMat = zeros(numOdor*numRecp,length(files));
        allSign = zeros(numOdor*numRecp,length(files));
        allfmin = zeros(length(files),1);
        
        for j0 = 1:length(files)
            temp = load(char(fullfile(folder,filesep,targetNames{i0},filesep,files{j0})));
%             allMat(:,j0)  = temp.wmin(:,1);
%             allSign(:,j0) = temp.wmin(:,2);

            allMat(:,j0)  = temp.wmin;
            allSign(:,j0) = temp.Sign;

            allfmin(j0) = temp.fmin;
        end
    end
    
      save([folder,filesep,'GcmiR0_N',num2str(numOdor),'_R',num2str(numRecp),'_S',num2str(sp),'_sig',num2str(sig),...
          '_basal',num2str(alpha),'_',dateInfo,'.mat'],'allMat','allfmin','allSign')
end

%% Extract and emerage data, reconstruction
close all
clear

% folder contain all th subfolders, default is the current folder
folder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding/recons_H50_N20_R9_S2_sig2_noise0.01_2018-09-26';
sFoler = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';  %folder to save data
% allFNames = dir(folder);
% folderFlag = [allFNames.isdir];
% subFolders = allFNames(folderFlag);
% subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
% dirNames = {subFolders.name};

N = 20;
M = 9;
sp = 2;
sig = 2;
noiseSig = 0.01;   %this should be matched to the data
H = 50;           % hidden layer size
% thd =  exp(-sig);   % threshold to detect
thd =  1;   % threshold to detect


allFile = dir(fullfile(folder,filesep,'*.mat'));
files = {allFile.name}';

allSp = 0.05:0.05:0.95;
% allSp = 0.1:0.1:1;



% for the old data
%{
allCrossEntr = [];
for i0 = 1:length(files)
    load(char(fullfile(folder,filesep,files{i0}))); 
    temp = zeros(size(summData,1),size(summData,2));
    for j0 = 1:size(summData,1)
        for k0 = 1:size(summData,2)
            temp(j0,k0) = summData(j0,k0).allTpr;
        end
    end
    allCrossEntr = [allCrossEntr,temp];
end
%}

allMSE = [];
for i0 = 1:length(files)
    load(char(fullfile(folder,filesep,files{i0}))); 
    temp = zeros(size(tr_all,1),size(tr_all,2));
    for j0 = 1:size(tr_all,1)
        for k0 = 1:size(tr_all,2)
            temp(j0,k0) = tr_all{j0,k0}(1,2);
        end
    end
%     allMSE = [allMSE,temp];
    allMSE(:,:,i0) = temp;
end


prfix = ['recons_N',num2str(N),'M',num2str(M),'sp',num2str(sp),'_noise',...
    num2str(noiseSig),'H',num2str(H),'thd',num2str(thd),'_',date,'.mat'];
save(fullfile(sFoler,prfix),'allMSE','allSp')


%% Extract reconstruction using tensorflow
dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding/Recon_inhi_N50M10_L6';
sFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

allFile = dir(fullfile(dFolder,filesep,'*.mat'));
files = {allFile.name}';

N = 50;
M = 10;
spar = 2;
sig = 2;
H = 100;   %hidden layer size
noise = 0.01;
repeat  = 40;
L = 6;   % number of hidden layer

% allSp = 0.05:0.05:1;   % the reversed order of sparsity of W
allSp = [0.05, 0.08, 0.1:0.02:0.18,0.25,0.3:0.1:0.7,0.82,0.84,0.88,0.9,0.95];     % the spontaneous activity


allTestLoss = nan(length(allSp),repeat);
allTrainLoss = nan(length(allSp),repeat);

for i0 = 1:length(files)
    st = '(?<= *_sp)[\d.]+(?=_)';
%     st = '(?<= *sp)[\d.]+(?=\.mat)';
    sp = str2num(char(regexp(files{i0},st,'match')));
    inx = find(round(100*allSp)==round(100*sp));
    
    load(char(fullfile(dFolder,filesep,files{i0})));
    allTrainLoss(inx,1:length(trainLost)) = trainLost;
    allTestLoss(inx,1:length(trainLost)) = testLost;
end

% save the summarized data
% dName = ['recons_inhi_N', num2str(N),'M', num2str(M),'sp',num2str(spar),...
%     'sig',num2str(sig),'ns',num2str(noise),'_loss_',date,'.mat'];
dName = ['recons_inhi_N', num2str(N),'M', num2str(M),'sp',num2str(spar),...
    'sig',num2str(sig),'ns',num2str(noise),'_L', num2str(L),'_loss_',date,'.mat'];
save(fullfile(sFolder,dName),'allTestLoss','allTrainLoss','N','M','noise','spar','sig','allSp','H')
% plot((1-allSp)',nanmean(allTestLoss,2))

% errorbar((1-allSp)',nanmean(allTestLoss,2), nanstd(allTestLoss,0,2))
errorbar(allSp',nanmean(allTestLoss,2), nanstd(allTestLoss,0,2))
xlabel('$r_0$','Interpreter','latex')
ylabel('error')


%% Extract and emerage data, classification
close all
clear

% folder contain all th subfolders, default is the current folder
folder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding/classify_N50_R10_S3_noiseSig0.01_nType3_2018-07-18';
sFoler = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';  %folder to save data
% allFNames = dir(folder);
% folderFlag = [allFNames.isdir];
% subFolders = allFNames(folderFlag);
% subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
% dirNames = {subFolders.name};

N = 50;
M = 10;      % number of receptors
sp = 3;      % sparsity of input
nType = 2;   % three different groups
noiseSig = 0.01;  % std of Gaussian noise


allFile = dir(fullfile(folder,filesep,'*.mat'));
files = {allFile.name}';

allSp = 0.05:0.05:1;

allMSE = [];
for i0 = 1:length(files)
    load(char(fullfile(folder,filesep,files{i0}))); 
    allMSE = [allMSE,summData.errorRate];
end

prfix = ['classify_N',num2str(N),'M',num2str(M),'sp',num2str(sp),'_noise',...
    num2str(noiseSig),'_nType',num2str(nType),'_',date,'.mat'];
save(fullfile(sFoler,prfix),'allCrossEntr','allSp')
% this program extract fmin and Wmin from the raw data, this happens when
% we run the simulation and forget to save the variables
%

close all
clear

%% load the data
% folder contain all th subfolders, default is the current folder
folder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/gcmi_sig/rawData';

allFNames = dir(folder);
folderFlag = [allFNames.isdir];
subFolders = allFNames(folderFlag);
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
dirNames = {subFolders.name};


% set the landmark in the folder names, so that we can just handle what we
% want. Right now we jut need the data folders that are from cluster

% sig =2;   %standard deviation of

% using regular expression to find the relavent folders
str_marker = '2018-07-22';   %folder with this string contains the data we need
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
        numOdor = str2num(char(regexp(targetNames{i0},s1,'match')));
        numRecp = str2num(char(regexp(targetNames{i0},s2,'match')));
        sp = str2num(char(regexp(targetNames{i0},s3,'match')));
        tp = char(regexp(targetNames{i0},s4,'match'));
        sig = str2num(char(regexp(targetNames{i0},s5,'match')));
        allMat = zeros(numOdor*numRecp,length(files));
        allfmin = zeros(length(files),1);
        
        for j0 = 1:length(files)
            load(char(fullfile(folder,filesep,targetNames{i0},filesep,files{j0})));
            
            % notice that the loaded file contain many irrelevant variabels
            allMat(:,j0)  = xmin;
            allfmin(j0) = fmin;
        end
    end
    
    
    
%     save(['u1m_N',num2str(numOdor),'_R',num2str(numRecp),...
%     '_S',num2str(S),'_sig',num2str(sig),'_',dateInfo,'.mat'],'allMat','allfmin')
%       save([folder,filesep,'spDepend_N',num2str(numOdor),'_R',num2str(numRecp),'_S',num2str(sp),'_sig',num2str(sig),...
%           '_',dateInfo,'.mat'],'allMat','allfmin')
      save([folder,filesep,tp,'int_N',num2str(numOdor),'_R',num2str(numRecp),'_S',num2str(sp),'_sig',num2str(sig),...
          '_',dateInfo,'.mat'],'allMat','allfmin')
end


% analysis of the optimal matrix
dFolder = '/Users/shan/Dropbox/olfactionProject/data/GcmiDifSigData';

% load the matrix
fName = 'N100_R30_S2_sig2_2018-05-18.mat';
load(fullfile(dFolder,fName));

% parameters
N = 100;
M = 30;
spar = 2;
sig =2;

repeats = 500;  %number of permuation

slect = 11;  %index to select one matrix, from 1 to 20
w = reshape(allMat(:,slect),[M,N]);



% ==============================================
% map the elements of original matrix into a stardard 
% Gaussian distribution and compare the eigenvalues
% ==============================================
[f,x] = ecdf(w(:));
f(end) = 1-1e-4;
[~, inx] = sort(w(:));   %index
newW = norminv(f(2:end));
w0 = zeros(size(w,1),size(w,2));
w0(inx) = newW;


corrType = 'Pearson';  %can be Spearman or Pearson
C0 = corr(w0,'type',corrType);
eigenOpt = eig(C0);

eigenPerm = zeros(N, repeats);  %store eigenvalues of permuated matrix
aveOpt = [mean(eigenOpt),std(eigenOpt)];

for i0 = 1:repeats
    newW = reshape(w0(randperm(M*N)),size(w0,1),size(w0,2));
    temp = corr(newW,'type',corrType);
    eigenPerm(:,i0) = eig(temp);
end
avePerm = [mean(eigenPerm,1);std(eigenPerm,0,1)];

largest = max(eigenPerm);  %largest eigen values


% plot
Bu = brewermap(11,'Blues');    %define blue colors
% compare the standard deviation of eigenvalues
figure
hold on
histogram(avePerm(2,:),'FaceColor',Bu(9,:),'Normalization','probability');
ah = gca;
ylim = ah.YLim;
plot([aveOpt(2); aveOpt(2)],ylim,'r--','LineWidth',2)
box on
legend('shuffle','original')
xlabel('$\sigma_{\lambda}$','Interpreter','latex')
ylabel('probability')

% compare the maxium eigenvalue
figure
hold on
histogram(largest,'FaceColor',Bu(9,:),'Normalization','probability');
ah = gca;
ylim = ah.YLim;
plot([max(eigenOpt); max(eigenOpt)],ylim,'r--','LineWidth',2)
box on
legend('shuffle','original')
xlabel('$\lambda_{max}$','Interpreter','latex')
ylabel('probability')


% ecdf of data and Gaussian
figure
plot(x,f)
xlabel('w')
ylabel('ecdf')

% standard Gaussian
pd = makedist('Normal');
X = -3:.1:3;
cdf_normal = cdf(pd,X);
figure
plot(X,cdf_normal)
xlabel('x')
ylabel('cdf')

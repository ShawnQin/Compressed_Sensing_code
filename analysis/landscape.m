clc
clear
lnw1 = [-4:0.1:4]';
lnw2 = lnw1;
lnc = normrnd(0,4,1,1e3);
lnc = 6*rand(1,1e3)-3;
m = [0.5 1 2];
figure('color','white')
for k = 1:3
    n = m(k);
r1 = (exp(lnw1)*exp(lnc)).^n./(1+(exp(lnw1)*exp(lnc)).^n);
r2 = (exp(lnw2)*exp(lnc)).^n./(1+(exp(lnw2)*exp(lnc)).^n);
for i = 1:length(lnw1)
    for j = 1:length(lnw2)
    r = [r1(i,:);r2(j,:)];
    p = kde(r,'rot');
    p1 = kde(r(1,:),'rot');
    p2 = kde(r(2,:),'rot');
    h(i,j) = entropy(p);
    h1(i) = entropy(p1);
    h2(j) = entropy(p2);
    hint(i,j) = h(i,j) - h1(i) - h2(j);
    end
end
subplot(2,3,k)
surf(lnw1,lnw2,h);
colormap jet;
xlabel('ln(weight1)');
ylabel('ln(weight2)');
zlabel('differential entropy')
set(gca,'FontSize',10)
axis tight;
title(strcat('1 odor 2 receptor -- n = ',32,num2str(m(k))));
subplot(2,3,k+3)
contourf(lnw1,lnw2,h,15);
xlabel('ln(weight1)');
ylabel('ln(weight2)');
set(gca,'FontSize',10)
end
saveas(gcf,'1_2_landscape.fig')
%%
figure('color','white')
m = [0.5];
for k = 1
    n = m(k);
lncc = [normrnd(0,1,1,1e3);normrnd(0,1,1,1e3)];
for i = 1:length(lnw1)
    for j = 1:length(lnw2)
        r = (exp([lnw1(i),lnw2(j)])*exp(lncc)).^n./((exp([lnw1(i),lnw2(j)])*exp(lncc)).^n+1);
        p = kde(r,'rot');
    h(i,j) = entropy(p);
    hint(i,j) = h(i,j) - h1(i)-h2(j);
    end
end
subplot(2,1,k)
surf(lnw1,lnw2,h);
xlabel('ln(weight1)');
ylabel('ln(weight2)');
zlabel('differential entropy')
set(gca,'FontSize',15,'FontWeight')
colormap jet;
axis tight
title(strcat('2 odor 1 receptor -- n=',32,num2str(m(k))));
subplot(2,1,k+1)
contourf(lnw1,lnw2,h,50);
xlabel('ln(weight1)');
ylabel('ln(weight2)');
set(gca,'FontSize',15)
end
saveas(gcf,'2_1_landscape.fig')
%%
clear
clc
lnw1 = [-5:0.2:5];
lnw2 = lnw1;
lnw3 = lnw1;
lnw4 = lnw1;
lncc = [normrnd(0,1,1,1e3);normrnd(0,4,1,1e3)];
for i = 1:length(lnw1);
    for j = 1:length(lnw2);
        for k = 1:length(lnw3);
            for m = 1:length(lnw4);
                r = exp([lnw1(i),lnw2(j);lnw3(k),lnw4(m)])*exp(lncc)./(1+exp([lnw1(i),lnw2(j);lnw3(k),lnw4(m)])*exp(lncc));
                col = HShannon_KDP_initialization(1);
                h(i,j,k,m) = HShannon_KDP_estimation(r,col);
            end
        end
    end
end
save('2_2_result_difSig2.mat')
%%
w1 = -3:0.05:3;
w2 = w1;
c = -8:0.05:8;
g = @(x) exp(2*x)./(1+exp(x)).^4;
for i = 1:length(w1)
    for j = 1:length(w1)
h(i,j) = 1/2*sum(0.05*normpdf(c,0,4).*log((g(w1(i)-c))+(g(w2(j)-c))));
    end
end

colormap jet




        

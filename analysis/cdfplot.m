subplot(1,2,1)

w = [-inf 0];
corrcoefMat = eye(2);
mu = [0 0];
sigma = [4 1];
covMat = corr2cov(sigma,corrcoefMat);
odorcon = exp(mvnrnd(mu,covMat,1e5));
resp = exp(w)*odorcon'./(1+exp(w)*odorcon');
x = 0:0.01:1;
y = x;
[a,b] = ecdf(resp);
plot(x,y,'LineWidth',2);
hold on
plot(b,a,'LineWidth',2);
hold on

w1 = -5:1:0;
for i = 1:6
    w = [w1(i) 0];
    resp = exp(w)*odorcon'./(1+exp(w)*odorcon');
    [a,b] = ecdf(resp);
    plot(b,a,'LineWidth',2);
    hold on
end

set(gca,'LineWidth',1,'FontSize',14);
xlabel('response');
ylabel('cdf');
title('cumulative distribution');

subplot(1,2,2)

x = -10:0.01:10;
y1 = normpdf(x,0,sigma(1));
plot(x,y1,'LineWidth',2);
hold on
plot(-Inf*ones(1,51),0:0.01:0.5,'--','LineWidth',2);
for i = 1:6
plot(w1(i)*ones(1,51),0:0.01:0.5,'--','LineWidth',2);
hold on
end
xlim([-10,10]);
set(gca,'LineWidth',1,'FontSize',14);


subplot(1,2,1)

sigma = 1:4;
for i = 1:4
    x = -5:0.01:5;
    y = normpdf(x,0,sigma(i));
    plot(x,y,'LineWidth',2);
    hold on
    xlim([-5,5]);
end
legend('sigma = 1','sigma = 2','sigma = 3','sigma = 4');
set(gca,'FontSize',14,'LineWidth',1);

    
subplot(1,2,2)

w = 0;
t = 0:0.01:1;
for i = 1:4
    odorcon = exp(normrnd(0,sigma(i),1,1e5));
    resp = exp(w)*odorcon./(1+exp(w)*odorcon);
    [a,b] = ecdf(resp);
    plot(b,a,'LineWidth',2);
    hold on
end
z = t;
plot(t,z,'LineWidth',1);   
legend('sigma = 1','sigma = 2','sigma = 3','sigma = 4');
set(gca,'FontSize',14,'LineWidth',1);


w = [0 -inf;-inf 0];
mu = [0 0];
sigma = [4 1];
covMat = corr2cov(sigma,corrcoefMat);
odorcon = exp(mvnrnd(mu,covMat,1e5));
resp = exp(w)*odorcon'./(1+exp(w)*odorcon');
scatter(resp(1),resp(2))
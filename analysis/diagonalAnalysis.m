sensmat = reshape(wmin,9,20);
[~,b] = sort(sensmat);
for i = 1:20
for j = 1:9
a(b(j,i),i) = j;
end
end
h=0;
for i = 1:9
[~,x] = sort(a(i,(h+1:end)));
a(:,[h+1:end]) = a(:,x(1:20-h)+h);
h1 = find(a(i,:)==1);
if isempty(h1) 
h1 = h;
end
h = max(h1);
end

allWmin = permute(allWmin,[1 3 2]);
wmin = reshape(allWmin(:,1,1),12,40);
x = 10.^([-8:0.5:2]);
for i = 1:9
testData(i,[21*(i-1)+1:21*i])=x;
end
resp = exp(wmin)*testData./(1+exp(wmin)*testData);
[coeff,score,latent] = pca(resp');
for i = 1:9
plot3(score(1+21*(i-1):21*i,1),score(1+21*(i-1):21*i,2),score(1+21*(i-1):21*i,3),'.-','LineWidth',1,'MarkerSize',15);
hold on
end
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
set(gca,'FontSize',20);

for i = 1:13
plot3(score(1+15*(i-1):15*i,1),score(1+15*(i-1):11*i,2),score(1+15*(i-1):15*i,3),'.-','LineWidth',1,'color',h(floor((i/2)^2+1),:));
hold on;scatter3(score(1+15*(i-1):15*i,1),score(1+11*(i-1):15*i,2),score(1+15*(i-1):15*i,3),[1:11].^2.*20,h(floor((i/2)^2+1),:),'.');hold on
end
 
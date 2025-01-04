clc
clear
load('2_2_result_difSig2.mat')
vidObj = VideoWriter('2_2_changing_w.avi');
vidObj.FrameRate = 2;
open(vidObj);
[x,y,z] = meshgrid(1:51,1:51,1:51);
x = reshape(x,1,length(lnw1).^3);
y = reshape(y,1,length(lnw1).^3);
z = reshape(z,1,length(lnw1).^3);
h1 = permute(h,[1 3 2 4]);
h2 = permute(h,[1 4 2 3]);
figure('color','white');
for i = 1:length(lnw1)
    for j = 1:length(lnw1)
    subplot(1,3,1)
    c = contourf(lnw1,lnw3,h1(:,:,i,j),15);
    clabel(c);
    caxis([-8,-1.5]);
    xlabel('concentration std 1 -- w21');
    ylabel('concentration std 1 -- w11')
    colormap jet
    colorbar
    title(strcat('lnw12 =',32,num2str(lnw1(i)),32,'lnw22 =',32,num2str(lnw1(j))));
    subplot(1,3,2)
    c1 = contourf(lnw1,lnw2,h(:,:,i,j),15);
    clabel(c1);
    caxis([-8,-1.5]);
    xlabel('concentration std 4 -- w12');
    ylabel('concentration std 1 -- w11')
    colormap jet
    colorbar
    title(strcat('lnw21 =',32,num2str(lnw1(i)),32,'lnw22 =',32,num2str(lnw1(j))));
    subplot(1,3,3)
    c2 = contourf(lnw1,lnw4,h2(:,:,i,j),15);
    clabel(c2);
    caxis([-8,-1.5]);
    xlabel('concentration std 4 -- w22');
    ylabel('concentration std 1 -- w11')
    colormap jet
    colorbar
    title(strcat('lnw12 =',32,num2str(lnw1(i)),32,'lnw21 =',32,num2str(lnw1(j))));
    set(gcf,'unit','normalized','position',[0.1,0.1,0.96,0.32]);
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    end
end
close(vidObj)
    
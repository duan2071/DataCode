function [undis_pic] = Undistortion(pic)

xy_parameter = load('xy_correction.txt');
z_parameter = load('z_correction.txt');
pic(:,1:2) = round((pic(:,1:2)+0.01) * 100);
% size(pic) (n, 3)
% figure;
% scatter3(pic(:,1),pic(:,2),pic(:,3),1);
% k = -1/200000;
k = -1/280000;
undis_pic = zeros(1,3);
num = 1;
for l1 = 1:240
    y = l1 - 120;
    for l2 = 1:320
        x = l2 - 160;
        x1 = round(x*(1 + k*x^2 + k*y^2));
        y1 = round(y*(1 + k*x^2 + k*y^2));
        y1 = y1 + 120;
        x1 = x1 + 160;
        for i = 1:size(pic,1)
            if pic(i,1) == x1 && pic(i,2) == y1
                undis_pic(num,:) = [l2,l1,pic(i,3)];
                num = num + 1;
            end
        end
    end
end
figure;
% size(pic)
% size(undis_pic)
scatter3(undis_pic(:,1),undis_pic(:,2),undis_pic(:,3),1);
for i = 1:size(undis_pic,1)
    undis_pic(i,1) = (undis_pic(i,3) * xy_parameter(1) + xy_parameter(2)) * (undis_pic(i,1) - 160.5) + 1.605;
    undis_pic(i,2) = (undis_pic(i,3) * xy_parameter(1) + xy_parameter(2)) * (undis_pic(i,2) - 120.5) + 1.205;
    undis_pic(i,3) = undis_pic(i,3) * z_parameter(1) + z_parameter(2);
end
p1 = pointCloud(undis_pic);
figure;
pcshow(p1);
axis equal
% set(gca,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 1 1]);
set(gca,'XDir','reverse');
set(gca,'Color','w');
axis off
title('图像去畸变、映射到世界坐标系（xyz）');
end


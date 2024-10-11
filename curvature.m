function [ans] = curvature(data)


%vec储存法向量
vec = zeros(size(data));
%q储存曲率
q = zeros(length(data),1);
 
k = 50;
% 搜索每个点的最邻近点
neighbors = knnsearch(data(:,1:3),data(:,1:3), 'k', k+1);
for i = 1:length(data)
    curtemp = neighbors(i,2:end);
    indpoint = data(curtemp,:);
    % 计算协方差并提取特征
    [v, c] = eig(cov(indpoint));
    %特征值按照升序排列1<2<3
    c = diag(c)';
    %计算特征值的总和
    z = sum(c);
    %计算曲率，用最小特征值除/特征值总和，这也是特征归一化
    p1 = c(:,1)/z;
    q(i,:) = abs(p1);
    
    %最小特征值对应的列向量就是法向量，dot是交叉相乘
    vec(i,:) = v(:,1)';
end

% min_val = min(q); % 获取最小值和最大值
% max_val = max(q);
% q = (q - min_val) / (max_val - min_val); % 计算MinMax归一化结果
q = 10*q;
ans = q;
 
% 读取x
x = data(1:5:end,1);
% 读取y
y = data(1:5:end,2);
% 读取z
z = data(1:5:end,3);
% uvw为法向量的三列
u = vec(1:5:end,1);
v = vec(1:5:end,2);
w = vec(1:5:end,3);
ptCloud = pointCloud(data);
figure;
pcshow(ptCloud.Location,q);
axis equal
% set(gca,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 1 1]);
set(gca,'XDir','reverse');
set(gca,'Color','w');
axis off
ch = colorbar('position',[0.78 0.3 0.01 0.4]);
ch.Title.String = 'curvature';
title("Curvature");

% hold on
% % 显示法向量
% quiver3(x,y,z,u,v,w);
% hold off

end


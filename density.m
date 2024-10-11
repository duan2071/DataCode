function [ans] = density (data)

% 读取点云数据
ptCloud = pointCloud(data);

% 设置邻域半径和最小邻居数
radius = 3; % 邻域半径
minNeighbors = 0; % 最小邻居数

% 计算点云点密度
densityy = zeros(ptCloud.Count, 1); % 存储点密度
for i = 1:ptCloud.Count
    % 获取当前点的坐标
    currentPoint = ptCloud.Location(i, :);
    
    % 在邻域半径内查找点
    [indices, ~] = findNeighborsInRadius(ptCloud, currentPoint, radius);
    
    % 判断邻居点数量是否满足要求
    if numel(indices) >= minNeighbors
        densityy(i) = numel(indices);
    end
end
densityy = densityy/10;
% 可视化点密度
% figure;
% pcshow(ptCloud.Location, densityy);
% ch = colorbar('position',[0.78 0.3 0.01 0.4]);
% title('density');
% axis equal
% % set(gca,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 1 1]);
% set(gca,'XDir','reverse');
% axis off
pointShow(data, "density", densityy);

ans = densityy;

end


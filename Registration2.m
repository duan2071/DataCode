% ICP algorithm
function [data_source, data_target, registration_matrix] = Registration2(undis_pic1, undis_pic2, selected_pic1, selected_pic2,mode)

% matrix = [0.755784808943431,-0.0205896197989088,-0.654496287328574,36.7301593095850; ...
%         -0.0181484460226928,0.998462961522678,-0.0523674362016957,6.48964059498661; ...
%         0.654568526952921,0.0514566033071924,0.754249469007952,-98.3134717536925; ...
%         0,0,0,1];
% matrix = [0.74226026050554	-0.0129564354557032	0.669986444978195	44.3078159129476;...
%         0.00469832619165178	0.999889113304315	0.014131058930218	-2.069300899814;...
%         -0.670095240548114	-0.00734110862027079	0.742238827277978	112.669463116907;...
%         0    0	0	1];
matrix = [0.726701217000399	0.0227979905979403	-0.686575263780191	37.9155705758904;...
    -0.00606306348287976	0.999623086219212	0.0267754506740173	-2.10705327358666;...
    0.686926910574446	-0.015295003180469	0.726565538961461	-101.79309071404;...
    0	0	0	1];

anti_matrix = zeros(4,4);
anti_matrix(1:3,1:3) = inv(matrix(1:3,1:3));
anti_matrix(1:3,4) = -anti_matrix(1:3,1:3)*matrix(1:3,4);
anti_matrix(4,4) = 1;

if isempty(selected_pic1) && isempty(selected_pic2)
    data_source = undis_pic1' * 100;
    data_target = undis_pic2' * 100;
else
    data_source = selected_pic1' * 100;
    data_target = selected_pic2' * 100;
end

%% 参数设置 ICP
kd = 1;
inlier_ratio = 0.9;
Tolerance = 0.0001;
step_Tolerance = 0.00001;
max_iteration = 20000;
show = 0;

%% 弿始ICP
if mode == 0
    T_final=eye(4,4);   %旋转矩阵初始倿
elseif mode == 1
    T_final = anti_matrix;
else
    T_final = matrix;
end
iteration=0;
Rf=T_final(1:3,1:3);
Tf=T_final(1:3,4);
data_source_old = data_source;
data_source=Rf*data_source+Tf*ones(1,size(data_source,2));    %初次更新点集（代表粗配准结果＿
err=1;
%% 迭代优化
while(1)
    iteration=iteration+1;
    if kd == 1
        %利用Kd-tree找出对应点集
        kd_tree = KDTreeSearcher(data_target','BucketSize',10);
        [index, dist] = knnsearch(kd_tree, data_source');
        z = 0;
    else
        %利用欧式距离找出对应点集
        k=size(data_source,2);
        for i = 1:k
            data_q1(1,:) = data_target(1,:) - data_source(1,i);    % 两个点集中的点x坐标之差
            data_q1(2,:) = data_target(2,:) - data_source(2,i);    % 两个点集中的点y坐标之差
            data_q1(3,:) = data_target(3,:) - data_source(3,i);    % 两个点集中的点z坐标之差
            distance = sqrt(data_q1(1,:).^2 + data_q1(2,:).^2 + data_q1(3,:).^2);  % 欧氏距离
            [dist(i), index(i)] = min(distance);   % 找到距离朿小的那个炿
        end
    end
    
    disp(['误差err=',num2str(mean(dist))]);
    disp(['迭代次数ieration=',num2str(iteration)]);
    err_rec(iteration) = mean(dist);
    
    % 按距离排序，只取前面占比为inlierratio内的点以应对外点
    [~, idx] = sort(dist);
    inlier_num = round(size(data_source,2)*inlier_ratio);
    idx = idx(1:inlier_num);
    data_source_temp = data_source(:,idx);
    dist = dist(idx);
    index = index(idx);
    data_mid = data_target(:,index);
    
    % 去中心化后SVD分解求解旋转矩阵与平移向釿
    [R_new, t_new] = rigidTransform3D(data_source_temp', data_mid');
    
    % 计算累计的旋转矩阵与平移向量
    Rf = R_new * Rf;
    Tf = R_new * Tf + t_new;
    
    %     更新点集
    %     data_source=R_new*data_source+t_new*ones(1,size(data_source,2));
    data_source=Rf*data_source_old+Tf*ones(1,size(data_source_old,2));
    
    % 显示中间结果
    if show == 1
        h = figure(7);
        scatter3(data_source(1,:),data_source(2,:),data_source(3,:),1);
        hold on;
        scatter3(data_target(1,:),data_target(2,:),data_target(3,:),1);
        hold off;
        daspect([1 1 1]);
        pause(0.1);
        drawnow
    end
    
    if err < Tolerance
        disp('—???????????????????????????');
        disp('情况1：优化结果已经达到目标，结束优化');
        break
    end
    if iteration > 1 && err_rec(iteration-1) - err_rec(iteration) < step_Tolerance && err_rec(iteration-1) - err_rec(iteration) > 0
        disp('—???????????????????????????');
        disp('情况2：迭代每丿步带来的优化到极限忼，结束优化');
        break
    end
    if iteration>=max_iteration
        disp('—???????????????????????????');
        disp('情况3：迭代已经达到最大次数，结束优化');
        break
    end
end

%% 计算朿后结果的误差
if kd == 1
    %利用Kd-tree找出对应点集
    kd_tree = KDTreeSearcher(data_target','BucketSize',10);
    [index, dist] = knnsearch(kd_tree, data_source');
else
    %利用欧式距离找出对应点集
    k=size(data_source,2);
    for i = 1:k
        data_q1(1,:) = data_target(1,:) - data_source(1,i);    % 两个点集中的点x坐标之差
        data_q1(2,:) = data_target(2,:) - data_source(2,i);    % 两个点集中的点y坐标之差
        data_q1(3,:) = data_target(3,:) - data_source(3,i);    % 两个点集中的点z坐标之差
        distance = sqrt(data_q1(1,:).^2 + data_q1(2,:).^2 + data_q1(3,:).^2);  % 欧氏距离
        [dist(i), index(i)] = min(distance);   % 找到距离朿小的那个炿
    end
end
disp(['朿终误差err=',num2str(mean(dist))]);
err_rec(iteration+1) = mean(dist);

%% 迭代优化过程中误差变化曲线
figure;
plot(0:iteration,err_rec);title('迭代误差变化曲线');
grid on

%% 朿后点云匹配的结果
figure;
scatter3(data_source(1,:),data_source(2,:),data_source(3,:),1);
hold on;
scatter3(data_target(1,:),data_target(2,:),data_target(3,:),1);
hold off;
daspect([1 1 1]);
set(gcf,'color','black');%黑色
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','w','ycolor','w','zcolor','w')
grid off

%disp('旋转矩阵的真值：');
%disp(T0);  %旋转矩阵真忿
disp('计算出的旋转矩阵＿');
T_final = [Rf,Tf];
T_final=[T_final;0,0,0,1];
disp(T_final);

registration_matrix = T_final;

if ~isempty(selected_pic1) && ~isempty(selected_pic2)
    data_source = undis_pic1' * 100;
    data_source = Rf*data_source+Tf*ones(1,size(data_source,2));
    data_target = undis_pic2' * 100;
end

end

%% 计算两个点集p，q的刚性变换参数，p和q的大小要丿臿
function [R, t] = rigidTransform3D(p, q)
n = cast(size(p, 1), 'like', p);
m = cast(size(q, 1), 'like', q);
% 去中心化
pmean = sum(p,1)/n;
p2 = bsxfun(@minus, p, pmean);
qmean = sum(q,1)/m;
q2 = bsxfun(@minus, q, qmean);
% 对协方差矩阵进行SVD分解
C = p2'*q2;
[U,~,V] = svd(C);
R = V*diag([1 1 sign(det(U*V'))])*U';
t = qmean' - R*pmean';
end

%% 对点云进行旋转与平移
function [data_q,T] = rotate(data,theta_x, theta_y, theta_z, t)
theta_x = theta_x/180*pi;
rot_x = [1 0 0;0 cos(theta_x) sin(theta_x);0 -sin(theta_x) cos(theta_x)];
theta_y = theta_y/180*pi;
rot_y = [cos(theta_y) 0 -sin(theta_y);0 1 0;sin(theta_y) 0 cos(theta_y)];
theta_z = theta_z/180*pi;
rot_z = [cos(theta_z) sin(theta_z) 0;-sin(theta_z) cos(theta_z) 0;0 0 1];

% 变换矩阵
T = rot_x*rot_y*rot_z;
T = [T,t'];
T=[T;0,0,0,1];
%化为齐次坐标
rows=size(data,2);
rows_one=ones(1,rows);
data=[data;rows_one];
%返回三维坐标
data_q=T*data;
data_q=data_q(1:3,:);
end


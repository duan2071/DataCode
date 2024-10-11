% palnar and curvatrue algorithm
function [data_source, data_target, registration_matrix] = Registration3(undis_pic1, undis_pic2, selected_pic1, selected_pic2,mode)

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

%% ��������
kd = 1;
max_iteration = 20000;
show = 0;
show_total_err = 1;
fine_registration = 1;
err_threshold_percent = 1.0 / 100;

%% ��ʼICP
if mode == 0
    T_final=eye(4,4);   %��ת�����ʼֵ
elseif mode == 1
    T_final = anti_matrix;
else
    T_final = matrix;
end

Rf=T_final(1:3,1:3);
Tf=T_final(1:3,4);

iteration = 0;
curvature_filter_ratio = 0.2; %0.2;
density_filter_ratio = 0.8; %0.8;
inlier_ratio = 0.5; % 0.5;
Tolerance = 0.0001;
step_Tolerance = 0.000001;

characteristic_target = zeros(size(data_target,2),5);
characteristic_source = zeros(size(data_source,2),5);

characteristic_target(:, 2) = density(data_target');
% characteristic_target(:, 1) = curvature(data_target');
characteristic_target(:, 3:5) = data_target';

characteristic_source(:, 2) = density(data_source');
% characteristic_source(:, 1) = curvature(data_source');
characteristic_source(:, 3:5) = data_source';

[~, target_density_sort_index] = sort(characteristic_target(:,2), 'descend');
[~, source_density_sort_index] = sort(characteristic_source(:,2), 'descend');

characteristic_source_density_sorted = characteristic_source(source_density_sort_index, :);
characteristic_target_density_sorted = characteristic_target(target_density_sort_index, :);

characteristic_target_density_filter_num = round(size(characteristic_target_density_sorted,1) * density_filter_ratio);
characteristic_source_density_filter_num = round(size(characteristic_source_density_sorted,1) * density_filter_ratio);

characteristic_target_density_filtered = characteristic_target_density_sorted(1:characteristic_target_density_filter_num,:);
characteristic_source_density_filtered = characteristic_source_density_sorted(1:characteristic_source_density_filter_num,:);

characteristic_target_density_filtered(:, 1) = curvature(characteristic_target_density_filtered(:, 3:5));
characteristic_source_density_filtered(:, 1) = curvature(characteristic_source_density_filtered(:, 3:5));


[~, target_curvature_sort_index] = sort(characteristic_target_density_filtered(:,1).*characteristic_target_density_filtered(:,2), 'descend');
[~, source_curvature_sort_index] = sort(characteristic_source_density_filtered(:,1).*characteristic_source_density_filtered(:,2), 'descend');

characteristic_source_curvature_sorted = characteristic_source_density_filtered(source_curvature_sort_index, :);
characteristic_target_curvature_sorted = characteristic_target_density_filtered(target_curvature_sort_index, :);

characteristic_target_curvature_filter_num = round(size(characteristic_target_curvature_sorted,1) * curvature_filter_ratio);
% characteristic_source_curvature_filter_num = round(size(characteristic_source_curvature_sorted,1) * curvature_filter_ratio);

characteristic_target_curvature_filtered = characteristic_target_curvature_sorted(1:characteristic_target_curvature_filter_num,:);
% characteristic_source_curvature_filtered = characteristic_source_curvature_sorted(1:characteristic_source_curvature_filter_num,:);
characteristic_source_curvature_threshold = characteristic_target_curvature_filtered(characteristic_target_curvature_filter_num, 1)*characteristic_target_curvature_filtered(characteristic_target_curvature_filter_num, 2);
characteristic_source_curvature_filtered = characteristic_source_curvature_sorted(characteristic_source_curvature_sorted(:,1).*characteristic_source_curvature_sorted(:,2) > characteristic_source_curvature_threshold, :);

pointShow(characteristic_target_density_filtered(:, 3:5), "data\_target\_density\_filtered", characteristic_source_density_filtered(:, 2));
pointShow(characteristic_source_density_filtered(:, 3:5), "data\_source\_density\_filtered", characteristic_source_density_filtered(:, 2));
pointShow(characteristic_target_curvature_filtered(:, 3:5), "data\_target\_curvature\_filtered", characteristic_target_curvature_filtered(:, 1).*characteristic_target_curvature_filtered(:, 2));
pointShow(characteristic_source_curvature_filtered(:, 3:5), "data\_source\_curvature\_filtered", characteristic_source_curvature_filtered(:, 1).*characteristic_source_curvature_filtered(:, 2));

curvature_filtered_tree = KDTreeSearcher(characteristic_target_curvature_filtered(:, 1:2), 'BucketSize', 10);
[curvature_preregist_index, curvature_preregist_dist] = knnsearch(curvature_filtered_tree, characteristic_source_curvature_filtered(:, 1:2));

[~, curvature_preregist_dist_index_sorted] = sort(curvature_preregist_dist);
inline_num = round(size(curvature_preregist_dist,1) * inlier_ratio);
curvature_preregist_dist_index_sorted = curvature_preregist_dist_index_sorted(1:inline_num);
characteristic_source_curvature_filtered_nearesrt = characteristic_source_curvature_filtered(curvature_preregist_dist_index_sorted, :);
curvature_preregist_index = curvature_preregist_index(curvature_preregist_dist_index_sorted);
characteristic_target_curvature_filtered_nearesrt_mid = characteristic_target_curvature_filtered(curvature_preregist_index, :);

data_source_curvature_filtered_old = characteristic_source_curvature_filtered_nearesrt(:, 3:5)';
data_source_curvature_filtered=Rf*data_source_curvature_filtered_old+Tf*ones(1,size(data_source_curvature_filtered_old,2));    %���θ��µ㼯���������׼�����
data_target_curvature_filtered_mid = characteristic_target_curvature_filtered_nearesrt_mid(:, 3:5);
data_source_regist = Rf*data_source+Tf*ones(1,size(data_source,2));

while(1)
    iteration=iteration+1;
    if kd == 1
        curvature_filtered_tree_regist = KDTreeSearcher(data_target_curvature_filtered_mid, 'BucketSize', 10);
        [curvature_regist_index, curvature_regist_dist] = knnsearch(curvature_filtered_tree_regist, data_source_curvature_filtered');
    else
        error('not implemented');
    end
    
    [~, curvature_regist_dist_index_sorted] = sort(curvature_regist_dist);
    inline_num = round(size(curvature_regist_dist,1) * inlier_ratio);
    curvature_regist_dist_index_sorted = curvature_regist_dist_index_sorted(1:inline_num);
    curvature_regist_index = curvature_regist_index(curvature_regist_dist_index_sorted);
    
    
    disp(['��ƥ�����err=',num2str(mean(curvature_regist_dist))]);
    disp(['��������ieration=',num2str(iteration)]);
    err_rec(iteration) = mean(curvature_regist_dist);
    
    if show_total_err == 1
        tree__temp = KDTreeSearcher(data_target', 'BucketSize', 10);
        [~, dist_temp] = knnsearch(tree__temp, data_source_regist');
        disp(['�����err=',num2str(mean(dist_temp))]);
        err_rec(iteration) = mean(dist_temp);
    end
    
    
    data_source_curvature_filtered_temp = data_source_curvature_filtered(:, curvature_regist_dist_index_sorted);
    data_target_curvature_filtered_mid_temp = data_target_curvature_filtered_mid(curvature_regist_index, :);
    
    [R_new, t_new] = rigidTransform3D(data_source_curvature_filtered_temp', data_target_curvature_filtered_mid_temp);
    
    Rf = R_new * Rf;
    Tf = R_new * Tf + t_new;
    
    data_source_curvature_filtered = Rf*data_source_curvature_filtered_old+Tf*ones(1,size(data_source_curvature_filtered_old,2));
    data_source_regist = Rf*data_source+Tf*ones(1,size(data_source,2));
    
    if show == 1
        figure(20);
        scatter3(data_source_curvature_filtered(1,:),data_source_curvature_filtered(2,:),data_source_curvature_filtered(3,:),1);
        hold on;
        scatter3(data_target_curvature_filtered_mid(:,1),data_target_curvature_filtered_mid(:,2),data_target_curvature_filtered_mid(:,3),1);
        hold off;
        figure(21);
        scatter3(data_source_regist(1,:),data_source_regist(2,:),data_source_regist(3,:),1);
        hold on;
        scatter3(data_target(1,:),data_target(2,:),data_target(3,:),1);
        hold off;
        daspect([1 1 1]);
        pause(0.1);
        drawnow
    end
    
    if err_rec(iteration) < Tolerance
        disp('-------------------------------');
        disp('���1���Ż�����Ѿ��ﵽĿ�꣬�����Ż�');
        break;
    end
    if iteration > 1 && (abs(err_rec(iteration-1) - err_rec(iteration)) < step_Tolerance || (err_rec(iteration) - err_rec(iteration-1)) / err_rec(iteration - 1) > err_threshold_percent)
        disp('-------------------------------');
        disp('���2������ÿ���������Ż������ޣ������Ż�');
        break;
    end
    if iteration>=max_iteration
        disp('-------------------------------');
        disp('���3�������Ѿ��ﵽ�������������Ż�');
        break;
    end
end

data_source_old = data_source;
data_source=Rf*data_source+Tf*ones(1,size(data_source,2));

if fine_registration == 1
    iteration_fine = 0;
    inlier_ratio = 0.05;% inlier_ratio = 1;   %��������ͼ6��Դƥ����
    Tolerance = 0.0001;
    step_Tolerance = 0.00001;
    
    while(1)
        iteration_fine=iteration_fine+1;
        if kd == 1
            %����Kd-tree�ҳ���Ӧ�㼯
            kd_tree = KDTreeSearcher(data_target','BucketSize',10);
            [index, dist] = knnsearch(kd_tree, data_source');
        else
            %����ŷʽ�����ҳ���Ӧ�㼯
            k=size(data_source,2);
            for i = 1:k
                data_q1(1,:) = data_target(1,:) - data_source(1,i);    % �����㼯�еĵ�x����֮��
                data_q1(2,:) = data_target(2,:) - data_source(2,i);    % �����㼯�еĵ�y����֮��
                data_q1(3,:) = data_target(3,:) - data_source(3,i);    % �����㼯�еĵ�z����֮��
                distance = sqrt(data_q1(1,:).^2 + data_q1(2,:).^2 + data_q1(3,:).^2);  % ŷ�Ͼ���
                [dist(i), index(i)] = min(distance);   % �ҵ�����cС���Ǹ���
            end
        end
        
        disp(['���err=',num2str(mean(dist))]);
        disp(['��������ieration=',num2str(iteration_fine)]);
        err_rec(iteration_fine+iteration) = mean(dist);
        
        % ����������ֻȡǰ��ռ��Ϊinlierratio�ڵĵ���Ӧ�����
        [~, idx] = sort(dist);
        inlier_num = round(size(data_source,2)*inlier_ratio);
        idx = idx(1:inlier_num);
        data_source_temp = data_source(:,idx);
        dist = dist(idx);
        index = index(idx);
        data_mid = data_target(:,index);
        
        % ȥ���Ļ���SVD�ֽ������ת������ƽ�����Y
        [R_new, t_new] = rigidTransform3D(data_source_temp', data_mid');
        
        % �����ۼƵ���ת������ƽ������
        Rf = R_new * Rf;
        Tf = R_new * Tf + t_new;
        
        %     ���µ㼯
        %     data_source=R_new*data_source+t_new*ones(1,size(data_source,2));
        data_source=Rf*data_source_old+Tf*ones(1,size(data_source_old,2));
        
        % ��ʾ�м���
        if show == 1
            figure(21);
            scatter3(data_source(1,:),data_source(2,:),data_source(3,:),1);
            hold on;
            scatter3(data_target(1,:),data_target(2,:),data_target(3,:),1);
            hold off;
            daspect([1 1 1]);
            pause(0.1);
            drawnow
        end
        
        if err_rec(iteration_fine+iteration) < Tolerance
            disp('-------------------------------');
            disp('���1���Ż�����Ѿ��ﵽĿ�꣬�����Ż�');
            break;
        end
        if iteration_fine > 1 && (abs(err_rec(iteration_fine-1+iteration) - err_rec(iteration_fine+iteration)) < step_Tolerance || (err_rec(iteration_fine+iteration) - err_rec(iteration_fine-1+iteration)) /  err_rec(iteration_fine-1+iteration) > err_threshold_percent)
            disp('-------------------------------');
            disp('���2������ÿ���������Ż������ޣ������Ż�');
            break;
        end
        if iteration_fine>=max_iteration
            disp('-------------------------------');
            disp('���3�������Ѿ��ﵽ�������������Ż�');
            break;
        end
    end
end
%% ����c���������
if kd == 1
    %����Kd-tree�ҳ���Ӧ�㼯
    kd_tree = KDTreeSearcher(data_target','BucketSize',10);
    [~, dist] = knnsearch(kd_tree, data_source');
else
    %����ŷʽ�����ҳ���Ӧ�㼯
    k=size(data_source,2);
    for i = 1:k
        data_q1(1,:) = data_target(1,:) - data_source(1,i);    % �����㼯�еĵ�x����֮��
        data_q1(2,:) = data_target(2,:) - data_source(2,i);    % �����㼯�еĵ�y����֮��
        data_q1(3,:) = data_target(3,:) - data_source(3,i);    % �����㼯�еĵ�z����֮��
        distance = sqrt(data_q1(1,:).^2 + data_q1(2,:).^2 + data_q1(3,:).^2);  % ŷ�Ͼ���
        [dist(i), index(i)] = min(distance);   % �ҵ�����cС���Ǹ���
    end
end

disp(['�c�����err=',num2str(mean(dist))]);
err_rec(iteration+1+iteration_fine) = mean(dist);

%% �����Ż����������仯����
figure;
plot(0:iteration+iteration_fine,err_rec);title('�������仯����');
grid on

%% �c�����ƥ��Ľ��
figure;
scatter3(data_source(1,:),data_source(2,:),data_source(3,:),1);
hold on;
scatter3(data_target(1,:),data_target(2,:),data_target(3,:),1);
hold off;
daspect([1 1 1]);
set(gcf,'color','black');%��ɫ
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','w','ycolor','w','zcolor','w')
grid off

%disp('��ת�������ֵ��');
%disp(T0);  %��ת�������
disp('���������ת�����');
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


%% ���������㼯p��q�ĸ��Ա任������p��q�Ĵ�СҪد�a
function [R, t] = rigidTransform3D(p, q)
n = cast(size(p, 1), 'like', p);
m = cast(size(q, 1), 'like', q);
% ȥ���Ļ�
pmean = sum(p,1)/n;
p2 = bsxfun(@minus, p, pmean);
qmean = sum(q,1)/m;
q2 = bsxfun(@minus, q, qmean);
% ��Э����������SVD�ֽ�
C = p2'*q2;
[U,~,V] = svd(C);
R = V*diag([1 1 sign(det(U*V'))])*U';
t = qmean' - R*pmean';
end
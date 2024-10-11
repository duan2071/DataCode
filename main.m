clear;
close all;

%% ����ͼ��
[pic1, pic2] = Loadpic();%pic1��pic2��ʽ?? n*3 ����(����δ��??
%% ͼ��ȥ����???ӳ�䵽��������ϵ��xyz??
undis_pic1 = Undistortion(pic1); %����߶ȱ任��mΪ��??
undis_pic2 = Undistortion(pic2);
%% ��׼����ɸ???�����??Ҫ��
selected_pic1 = Select(undis_pic1,[],[],[]);
selected_pic2 = Select(undis_pic2,[],[],[]);

% load 20240607_selected.mat

%% ͼ����׼
tic;
[data_source1, data_target1, registration_matrix1] = Registration3(undis_pic1, undis_pic2, selected_pic1, selected_pic2,0);
toc

data_target4 = (registration_matrix1(1:3,1:3))\(data_target1 - registration_matrix1(1:3,4)*ones(1,size(data_target1,2)));


%% չʾͼ��
data_source_original = undis_pic1 * 100;
figure;
scatter3(data_source_original(:,1),data_source_original(:,2),data_source_original(:,3),1);
hold on;
scatter3(data_target4(1,:),data_target4(2,:),data_target4(3,:),1);
% hold on;
% scatter3(data_target5(1,:),data_target5(2,:),data_target5(3,:),1,data_target5(3,:));
% hold on;
% scatter3(data_target6(1,:),data_target6(2,:),data_target6(3,:),1);
hold off;
set(gcf,'color','k')
% set(gca,'XDir','reverse');
axis equal
axis off

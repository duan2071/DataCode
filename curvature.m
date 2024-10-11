function [ans] = curvature(data)


%vec���淨����
vec = zeros(size(data));
%q��������
q = zeros(length(data),1);
 
k = 50;
% ����ÿ��������ڽ���
neighbors = knnsearch(data(:,1:3),data(:,1:3), 'k', k+1);
for i = 1:length(data)
    curtemp = neighbors(i,2:end);
    indpoint = data(curtemp,:);
    % ����Э�����ȡ����
    [v, c] = eig(cov(indpoint));
    %����ֵ������������1<2<3
    c = diag(c)';
    %��������ֵ���ܺ�
    z = sum(c);
    %�������ʣ�����С����ֵ��/����ֵ�ܺͣ���Ҳ��������һ��
    p1 = c(:,1)/z;
    q(i,:) = abs(p1);
    
    %��С����ֵ��Ӧ�����������Ƿ�������dot�ǽ������
    vec(i,:) = v(:,1)';
end

% min_val = min(q); % ��ȡ��Сֵ�����ֵ
% max_val = max(q);
% q = (q - min_val) / (max_val - min_val); % ����MinMax��һ�����
q = 10*q;
ans = q;
 
% ��ȡx
x = data(1:5:end,1);
% ��ȡy
y = data(1:5:end,2);
% ��ȡz
z = data(1:5:end,3);
% uvwΪ������������
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
% % ��ʾ������
% quiver3(x,y,z,u,v,w);
% hold off

end


function [selected_pic] = Select(pic,x_range,y_range,z_range)

if isempty(x_range)
    x_range = [min(pic(:,1)),max(pic(:,1))];
end
if isempty(y_range)
    y_range = [min(pic(:,2)),max(pic(:,2))];
end
if isempty(z_range)
    z_range = [min(pic(:,3)),max(pic(:,3))];
end

num = 1;
for i = 1:size(pic,1)   %返回pic的行数
    if pic(i,1) >= x_range(1) && pic(i,1) <= x_range(2) && ...
            pic(i,2) >= y_range(1) && pic(i,2) <= y_range(2) && ...
            pic(i,3) >= z_range(1) && pic(i,3) <= z_range(2)
        selected_pic(num,:) = pic(i,:);
        num = num + 1;
    end
end
           
end


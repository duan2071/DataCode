function []=pointShow(data, title_str, cmap, position)
    if(~exist('position', 'var'))
        position = [0.78 0.3 0.01 0.4];
    end
    % ptCloud = pointCloud(data);
    figure;
    if(exist('cmap', 'var'))
        pcshow(data,cmap);
    else
        pcshow(data);
    end
    axis equal
    % set(gca,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 1 1]);
    set(gca,'XDir','reverse');
    set(gca,'Color','w');
    axis off
    ch = colorbar('position', position);
    ch.Title.String = title_str;
    title(title_str);
end
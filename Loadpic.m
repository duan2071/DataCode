function [pic1, pic2] = Loadpic()

[fileName,pathName]=uigetfile("*.txt","�ı��ļ�");
pic1=importdata([pathName,fileName]);
pic1=pic1(:,1:3);
% figure;
% scatter3(pic1(:,1),pic1(:,2),pic1(:,3),1);
[fileName,pathName]=uigetfile("*.txt","�ı��ļ�");
pic2=importdata([pathName,fileName]);
pic2=pic2(:,1:3);
% figure;
% scatter3(pic2(:,1),pic2(:,2),pic2(:,3),1);

end


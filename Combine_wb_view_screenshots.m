clear all; close all; clc

D = dir('Results/wb_view_screenshots/Group*');
I = imread(sprintf('Results/wb_view_screenshots/%s',D(1).name));
I = imcrop(I,[0 0 size(I,2)/2 size(I,1)]);
for n = 2:length(D)
    I2 = imread(sprintf('Results/wb_view_screenshots/%s',D(n).name));
    I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
    I = cat(2,I,I2);
end
imshow(I);
print(gcf,'Results/Group_maps_screenshots','-dpng','-r300');

D = dir('Results/wb_view_screenshots/SNR_Group*');
I = imread(sprintf('Results/wb_view_screenshots/%s',D(1).name));
I = imcrop(I,[0 0 size(I,2)/2 size(I,1)]);
for n = 2:length(D)
    I2 = imread(sprintf('Results/wb_view_screenshots/%s',D(n).name));
    I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
    I = cat(2,I,I2);
end
imshow(I);
print(gcf,'Results/Group_maps_matched_SNR_screenshots','-dpng','-r300');

D = dir('Results/wb_view_screenshots/Subject_Green_*');
I = imread(sprintf('Results/wb_view_screenshots/%s',D(1).name));
I = imcrop(I,[0 0 size(I,2)/2 size(I,1)]);
for n = 2:length(D)
    I2 = imread(sprintf('Results/wb_view_screenshots/%s',D(n).name));
    I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
    I = cat(2,I,I2);
end
imshow(I);
print(gcf,'Results/Subject_green_Individual_screenshots','-dpng','-r300');

D = dir('Results/wb_view_screenshots/Subject_Red_*');
I = imread(sprintf('Results/wb_view_screenshots/%s',D(1).name));
I = imcrop(I,[0 0 size(I,2)/2 size(I,1)]);
for n = 2:length(D)
    I2 = imread(sprintf('Results/wb_view_screenshots/%s',D(n).name));
    I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
    I = cat(2,I,I2);
end
imshow(I);
print(gcf,'Results/Subject_red_Individual_screenshots','-dpng','-r300');

D = dir('Results/wb_view_screenshots/Subject_from_group_Green_*');
I = imread(sprintf('Results/wb_view_screenshots/%s',D(1).name));
I = imcrop(I,[0 0 size(I,2)/2 size(I,1)]);
for n = 2:length(D)
    I2 = imread(sprintf('Results/wb_view_screenshots/%s',D(n).name));
    I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
    I = cat(2,I,I2);
end
imshow(I);
print(gcf,'Results/Subject_green_Group_screenshots','-dpng','-r300');

D = dir('Results/wb_view_screenshots/Subject_from_group_Red_*');
I = imread(sprintf('Results/wb_view_screenshots/%s',D(1).name));
I = imcrop(I,[0 0 size(I,2)/2 size(I,1)]);
for n = 2:length(D)
    I2 = imread(sprintf('Results/wb_view_screenshots/%s',D(n).name));
    I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
    I = cat(2,I,I2);
end
imshow(I);
print(gcf,'Results/Subject_red_Group_screenshots','-dpng','-r300');

D = dir('Results/wb_view_screenshots/Subject_Outlier_*');
I = imread(sprintf('Results/wb_view_screenshots/%s',D(1).name));
I = imcrop(I,[0 0 size(I,2)/2 size(I,1)]);
for n = 2:length(D)
    I2 = imread(sprintf('Results/wb_view_screenshots/%s',D(n).name));
    I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
    I = cat(2,I,I2);
end
imshow(I);
print(gcf,'Results/Subject_outlier_Individual_screenshots','-dpng','-r300');

D = dir('Results/wb_view_screenshots/Subject_from_group_Outlier_*');
I = imread(sprintf('Results/wb_view_screenshots/%s',D(1).name));
I = imcrop(I,[0 0 size(I,2)/2 size(I,1)]);
for n = 2:length(D)
    I2 = imread(sprintf('Results/wb_view_screenshots/%s',D(n).name));
    I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
    I = cat(2,I,I2);
end
imshow(I);
print(gcf,'Results/Subject_outlier_Group_screenshots','-dpng','-r300');

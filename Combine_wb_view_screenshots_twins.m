clear all; close all; clc

subs = [105923 114823 125525 130518 137128 144226 146129 158035 169343 177746 185442 192439 195041 200614 204521 562345 601127 783462 859671 861456];
Twins = [3 12; 4 16; 5 19; 7 18; 8 10; 11 17; 13 15; 14 20];
NotTwins = setdiff(1:20,Twins(:));

% Create Find-The-Biggest figure
Ia = [];
for t = 1:4
    I1 = imread(sprintf('Results/wb_view_screenshots/FTB_sub%02d.png',Twins(t,1)));
    I2 = imread(sprintf('Results/wb_view_screenshots/FTB_sub%02d.png',Twins(t,2)));
    I1 = cat(1,I1,I2); Ia = cat(2,Ia,I1);
end
fh = figure; imshow(Ia,'border','tight')
rectangle('Position',[1 1 2250-1 1658*2-1],'LineWidth',2);
rectangle('Position',[2250+1 1 2250-1 1658*2-1],'LineWidth',2);
rectangle('Position',[2*2250+1 1 2250-1 1658*2-1],'LineWidth',2);
rectangle('Position',[3*2250+1 1 2250-1 1658*2-1],'LineWidth',2);
Ia = getframe(fh); 

Ib = [];
for t = 5:8
    I1 = imread(sprintf('Results/wb_view_screenshots/FTB_sub%02d.png',Twins(t,1)));
    I2 = imread(sprintf('Results/wb_view_screenshots/FTB_sub%02d.png',Twins(t,2)));
    I1 = cat(1,I1,I2); Ib = cat(2,Ib,I1);
end
fh = figure; imshow(Ib,'border','tight')
rectangle('Position',[1 1 2250-1 1658*2-1],'LineWidth',2);
rectangle('Position',[2250+1 1 2250-1 1658*2-1],'LineWidth',2);
rectangle('Position',[2*2250+1 1 2250-1 1658*2-1],'LineWidth',2);
rectangle('Position',[3*2250+1 1 2250-1 1658*2-1],'LineWidth',2);
Ib = getframe(fh); 

Ic = [];
for t = 1:4
    I1 = imread(sprintf('Results/wb_view_screenshots/FTB_sub%02d.png',NotTwins(t)));
    Ic = cat(2,Ic,I1);
end
fh = figure; imshow(Ic,'border','tight')
rectangle('Position',[1 1 2250-1 1658-1],'LineWidth',2);
rectangle('Position',[2250+1 1 2250-1 1658-1],'LineWidth',2);
rectangle('Position',[2*2250+1 1 2250-1 1658-1],'LineWidth',2);
rectangle('Position',[3*2250+1 1 2250-1 1658-1],'LineWidth',2);
Ic = getframe(fh); 

I = cat(1,Ia.cdata,Ib.cdata,Ic.cdata);
imshow(I,'border','tight');
print(gcf,'Results/FTB_screenshots','-dpng','-r300');

% Create Overlap figure
Ia = [];
for t = 1:4
    I1 = imread(sprintf('Results/wb_view_screenshots/Overlap_sub%02d.png',Twins(t,1)));
    I2 = imread(sprintf('Results/wb_view_screenshots/Overlap_sub%02d.png',Twins(t,2)));
    I1 = cat(1,I1,I2); Ia = cat(2,Ia,I1);
end
fh = figure; imshow(Ia,'border','tight')
rectangle('Position',[1 1 2250-1 1658*2-1],'LineWidth',2);
rectangle('Position',[2250+1 1 2250-1 1658*2-1],'LineWidth',2);
rectangle('Position',[2*2250+1 1 2250-1 1658*2-1],'LineWidth',2);
rectangle('Position',[3*2250+1 1 2250-1 1658*2-1],'LineWidth',2);
Ia = getframe(fh); 

Ib = [];
for t = 5:8
    I1 = imread(sprintf('Results/wb_view_screenshots/Overlap_sub%02d.png',Twins(t,1)));
    I2 = imread(sprintf('Results/wb_view_screenshots/Overlap_sub%02d.png',Twins(t,2)));
    I1 = cat(1,I1,I2); Ib = cat(2,Ib,I1);
end
fh = figure; imshow(Ib,'border','tight')
rectangle('Position',[1 1 2250-1 1658*2-1],'LineWidth',2);
rectangle('Position',[2250+1 1 2250-1 1658*2-1],'LineWidth',2);
rectangle('Position',[2*2250+1 1 2250-1 1658*2-1],'LineWidth',2);
rectangle('Position',[3*2250+1 1 2250-1 1658*2-1],'LineWidth',2);
Ib = getframe(fh); 

Ic = [];
for t = 1:4
    I1 = imread(sprintf('Results/wb_view_screenshots/Overlap_sub%02d.png',NotTwins(t)));
    Ic = cat(2,Ic,I1);
end
fh = figure; imshow(Ic,'border','tight')
rectangle('Position',[1 1 2250-1 1658-1],'LineWidth',2);
rectangle('Position',[2250+1 1 2250-1 1658-1],'LineWidth',2);
rectangle('Position',[2*2250+1 1 2250-1 1658-1],'LineWidth',2);
rectangle('Position',[3*2250+1 1 2250-1 1658-1],'LineWidth',2);
Ic = getframe(fh); 

I = cat(1,Ia.cdata,Ib.cdata,Ic.cdata);
imshow(I,'border','tight');
print(gcf,'Results/Overlap_screenshots','-dpng','-r300');

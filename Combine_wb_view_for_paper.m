clear all; close all; clc

% Figure of all group maps
Dgroup = dir('Results/wb_view_screenshots/Group*');
Names = {'L-FPN 1 (control 1)','ON 1 (visual 1)','D-FPN 1 (attention 1)','PN (somatomotor)','ON 2 (visual 2)','ON 3 (visual 3)',...
    'M-FPN 1 (default 1)','M-FPN 2 (default 2)','L-FPN 2 (control 2)','L-FPN 3 (control 3)','D-FPN 2 (attention 2)','L-FPN 4 (control 4)'};
t = tiledlayout(3,4,'TileSpacing','Compact','Padding','None');
nexttile
for n = 1:length(Dgroup)
    I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgroup(n).name));
    I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
    imshow(I2,'border','tight'); 
    title(Names{n},'FontSize',14);
    if n<12
        nexttile
    end
    
end
set(gcf,'Position',[1 1 661 687],'PaperPositionMode','auto')
print(gcf,'Results/Group_maps_for_paper','-dpng','-r300');

% Figure focused on maps 7 and 8
Dsnr = dir('Results/wb_view_screenshots/SNR_Group*');
Dgreen = dir('Results/wb_view_screenshots/Subject_Green_*');
Dred = dir('Results/wb_view_screenshots/Subject_Red_*');
Dgreengrp = dir('Results/wb_view_screenshots/Subject_from_group_Green_*');
Dredgrp = dir('Results/wb_view_screenshots/Subject_from_group_Red_*');

% Map 7
t = tiledlayout(3,2,'TileSpacing','Compact','Padding','None');
nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgroup(7).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title('Group','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dsnr(7).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Group','(matched SNR)'},'FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgreen(7).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 195041','(individual run)'},'color','g','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgreengrp(7).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 195041','(group run)'},'color','g','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dred(7).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 125525','(individual run)'},'color','r','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dredgrp(7).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 125525','(group run)'},'color','r','FontSize',14); 

set(gcf,'Position',[1 1 350 687],'PaperPositionMode','auto')
print(gcf,'Results/Comparison_maps_for_paper_A','-dpng','-r300');

% Map 8 
t = tiledlayout(3,2,'TileSpacing','Compact','Padding','None');
nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgroup(8).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title('Group','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dsnr(8).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Group','(matched SNR)'},'FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgreen(8).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 195041','(individual run)'},'color','g','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgreengrp(8).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 195041','(group run)'},'color','g','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dred(8).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 125525','(individual run)'},'color','r','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dredgrp(8).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 125525','(group run)'},'color','r','FontSize',14); 

set(gcf,'Position',[1 1 350 687],'PaperPositionMode','auto')
print(gcf,'Results/Comparison_maps_for_paper_B','-dpng','-r300');

% % Map 7&8 
% t = tiledlayout(1,2,'TileSpacing','Compact','Padding','None');
% nexttile
% imshow(t1.cdata,'border','tight'); title(sprintf('A. %s',Names{7}),'FontSize',14);
% rectangle('Position',[1 1 700 1374],'EdgeColor','k','LineWidth',1);
% nexttile
% imshow(t2.cdata,'border','tight'); title(sprintf('B. %s',Names{8}),'FontSize',14);
% rectangle('Position',[1 1 700 1374],'EdgeColor','k','LineWidth',1);
% set(gcf,'Position',[1 1 661 687],'PaperPositionMode','auto')
% 

% OHBM ABSTRACT Figure focused on maps 7 and 8
Dsnr = dir('Results/wb_view_screenshots/SNR_Group*');
Dgreen = dir('Results/wb_view_screenshots/Subject_Green_*');
Dred = dir('Results/wb_view_screenshots/Subject_Red_*');
Dgreengrp = dir('Results/wb_view_screenshots/Subject_from_group_Green_*');
Dredgrp = dir('Results/wb_view_screenshots/Subject_from_group_Red_*');

% Map 7
t = tiledlayout(2,6,'TileSpacing','Compact','Padding','None');
nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgroup(7).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title('Group','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dsnr(7).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Group','(matched SNR)'},'FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgreen(7).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 195041','(individual run)'},'color','g','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgreengrp(7).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 195041','(group run)'},'color','g','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dred(7).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 125525','(individual run)'},'color','r','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dredgrp(7).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 125525','(group run)'},'color','r','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgroup(8).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title('Group','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dsnr(8).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Group','(matched SNR)'},'FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgreen(8).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 195041','(individual run)'},'color','g','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dgreengrp(8).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 195041','(goup run)'},'color','g','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dred(8).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 125525','(individual run)'},'color','r','FontSize',14); nexttile

I2 = imread(sprintf('Results/wb_view_screenshots/%s',Dredgrp(8).name)); I2 = imcrop(I2,[0 0 size(I2,2)/2 size(I2,1)]);
imshow(I2,'border','tight'); title({'Subject 125525','(group run)'},'color','r','FontSize',14); 

set(gcf,'Position',[1 1 1050 475],'PaperPositionMode','auto')
print(gcf,'Results/Comparison_maps_for_abstract','-dpng','-r300');

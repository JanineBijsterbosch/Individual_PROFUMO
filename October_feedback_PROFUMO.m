clear all; close all; clc

addpath(genpath('~/Dropbox/Matlab/fieldtrip-20190819/'),'-END')
addpath('~/Dropbox/Matlab/','-END')
load('PFM_keep.mat')

%% [matt] Group PROFUMO at different dimensionalitiees
DimCheck = zeros(12,4);
grp20 = ft_read_cifti('PROFUMO/0All.pfm/Results.ppp/Maps/Group.dscalar.nii'); grp20 = dscalar2double(grp20,1); 
%grp12 = ft_read_cifti('PROFUMO/0All.pfm/Results.ppp/Maps/Group_D12.dscalar.nii'); grp12 = dscalar2double(grp12,1);
grp30 = ft_read_cifti('PROFUMO/0All.pfm/Results.ppp/Maps/Group_D30.dscalar.nii'); grp30 = dscalar2double(grp30,1);
grp40 = ft_read_cifti('PROFUMO/0All.pfm/Results.ppp/Maps/Group_D40.dscalar.nii'); grp40 = dscalar2double(grp40,1);
grp50 = ft_read_cifti('PROFUMO/0All.pfm/Results.ppp/Maps/Group_D50.dscalar.nii'); grp50 = dscalar2double(grp50,1);
   
%assign = munkres(1-corr(grp12,grp20)); [i,~] = find(assign); grp12 = grp12(:,i); grp12 = grp12(:,PFM_keep_group);
assign = munkres(1-corr(grp30,grp20)); [i,~] = find(assign); grp30 = grp30(:,i); grp30 = grp30(:,PFM_keep_group);
assign = munkres(1-corr(grp40,grp20)); [i,~] = find(assign); grp40 = grp40(:,i); grp40 = grp40(:,PFM_keep_group);
assign = munkres(1-corr(grp50,grp20)); [i,~] = find(assign); grp50 = grp50(:,i); grp50 = grp50(:,PFM_keep_group);
grp20 = grp20(:,PFM_keep_group);    
%A = corr(grp20,grp12); DimCheck(:,1) = A(eye(12)==1); clear A
A = corr(grp20,grp30); DimCheck(:,2) = A(eye(12)==1); clear A
A = corr(grp20,grp40); DimCheck(:,3) = A(eye(12)==1); clear A
A = corr(grp20,grp50); DimCheck(:,4) = A(eye(12)==1); clear A
figure;
% swarmchart(ones(1,size(DimCheck,1)),DimCheck(:,1),20,'cyan','filled'); hold on; 
% plot([0.75 1.25],[nanmean(DimCheck(:,1)) nanmean(DimCheck(:,1))],'k'); plot([0.75 1.25],[nanmedian(DimCheck(:,1)) nanmedian(DimCheck(:,1))],'k','linewidth',2);
swarmchart(ones(1,size(DimCheck,1))+0,DimCheck(:,2),20,'magenta','filled'); hold on; 
plot([0.75 1.25],[nanmean(DimCheck(:,2)) nanmean(DimCheck(:,2))],'k'); plot([0.75 1.25],[nanmedian(DimCheck(:,2)) nanmedian(DimCheck(:,2))],'k','linewidth',2);
swarmchart(ones(1,size(DimCheck,1))+1,DimCheck(:,3),20,'green','filled'); hold on;
plot([1.75 2.25],[nanmean(DimCheck(:,3)) nanmean(DimCheck(:,3))],'k'); plot([1.75 2.25],[nanmedian(DimCheck(:,3)) nanmedian(DimCheck(:,3))],'k','linewidth',2);
swarmchart(ones(1,size(DimCheck,1))+2,DimCheck(:,4),20,'blue','filled'); hold on; 
plot([2.75 3.25],[nanmean(DimCheck(:,4)) nanmean(DimCheck(:,4))],'k'); plot([2.75 3.25],[nanmedian(DimCheck(:,4)) nanmedian(DimCheck(:,4))],'k','linewidth',2);
set(gca,'xtick',[1 2 3],'xticklabel',{'D=30','D=40','D=50'},'FontSize',14);
title('Group map spatial correlations across dimensionalities','FontSize',14)
ylabel('Spatial correlation with D=20','FontSize',14)
print('Results/October_feedback_Dimensionalitiy','-dpng','-r300');

%% [rezvan] Would be good to look at distribution of spatial PROFUMO weights between subject and group maps (either as scatter plot or comparing two distributions)
subject_maps_separate_runs = ft_read_cifti('Results/Maps/subject_maps_separate_runs.dtseries.nii'); subject_maps_separate_runs = subject_maps_separate_runs.dtseries;
subject_maps_group_run = ft_read_cifti('Results/Maps/subject_maps_group_run.dtseries.nii'); subject_maps_group_run = subject_maps_group_run.dtseries;
s = [12 19]; figure; set(gcf,'position',[620 315 1181 732],'PaperPositionMode','auto')
for sI = 1:length(s)
    subplot(2,2,sI)
    scatter(subject_maps_separate_runs(:,(s(sI)-1)*length(PFM_keep_group)+1:s(sI)*length(PFM_keep_group)),subject_maps_group_run(:,(s(sI)-1)*length(PFM_keep_group)+1:s(sI)*length(PFM_keep_group)));
    axis([-6 16 -6 16])
    xlabel('spatial PROFUMO weights separate subject runs')
    ylabel('spatial PROFUMO weights subject map from group run')
    title(sprintf('subject %d',subs(s(sI))))
end
I1 = find(subject_maps_separate_runs>0.1 | subject_maps_separate_runs<-0.1);
I2 = find(subject_maps_group_run>0.1 | subject_maps_group_run<-0.1);
I = intersect(I1,I2);
subplot(2,2,3:4); 
histogram(subject_maps_group_run(I)); hold on; histogram(subject_maps_separate_runs(I));
ylabel('spatial PROFUMO weights')
legend({'subject maps from group run','separate subject runs'})
title('Distribution of PROFUMO spatial weights (exlcuding below abs(0.1) for visual clarity)')
print('Results/October_feedback_PROFUMO','-dpng','-r300');

%% [steve] Would be interesting to look at the most unique (outlier) subject to see if they still have good individual-group similarity and good correspondence
R1 = corrs(1,:); R1 = reshape(R1,12,20);
[~,outlier_mean] = min(nanmean(R2));
[~,outlier_median] = min(nanmedian(R2));
if outlier_mean == outlier_median
    outlier_sub = outlier_mean;
end
fprintf('Outlier subject: %d \n',subs(outlier_sub));
example = ft_read_cifti('OLD/Overlap.dtseries.dtseries.nii');
example.dtseries = subject_maps_separate_runs(:,(outlier_sub-1)*length(PFM_keep_group)+1:outlier_sub*length(PFM_keep_group)); example.hdr.dim(7) = length(PFM_keep_group); example.time = 1:length(PFM_keep_group);
ft_write_cifti(sprintf('Results/Maps/Example_Subject_OUTLIER_%d',subs(outlier_sub)),example,'parameter','dtseries');
example.dtseries = subject_maps_group_run(:,(outlier_sub-1)*length(PFM_keep_group)+1:outlier_sub*length(PFM_keep_group)); example.hdr.dim(7) = length(PFM_keep_group); example.time = 1:length(PFM_keep_group);
ft_write_cifti(sprintf('Results/Maps/Example_Subject_OUTLIER_%d_from_group',subs(outlier_sub)),example,'parameter','dtseries');



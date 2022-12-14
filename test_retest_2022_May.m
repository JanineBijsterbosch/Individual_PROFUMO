clear all; close all; clc

addpath(genpath('~/Dropbox/Matlab/fieldtrip-20190819/'),'-END')
addpath('~/Dropbox/Matlab/','-END')
addpath(genpath('~/Dropbox/Matlab/FSLNets'),'-END')
%addpath(genpath('~/Dropbox/Matlab/Violinplot-Matlab-master'))
addpath(genpath('~/Dropbox/Matlab/ICC'),'-END')
addpath('~/Dropbox/Matlab/hline_vline')

subs = [105923 114823 125525 130518 137128 144226 146129 158035 169343 177746 185442 192439 195041 200614 204521 562345 601127 783462 859671 861456];

% Test-retest reliability and group-subject similarity
trt_corrs = nan(size(subs,2),20); trt_labels = cell(size(trt_corrs));
trt_grpcorrs = nan(size(subs,2)*2,20); trt_grplabels = cell(size(trt_grpcorrs));
for s = 1:size(subs,2)
    fprintf('Running subject %d (%d)\n',s,subs(s))
    split1 = ft_read_cifti(sprintf('PROFUMO/2022_May/%d_split1.pfm/Results.ppp/Maps/Group.dscalar.nii',subs(s)));
    split2 = ft_read_cifti(sprintf('PROFUMO/2022_May/%d_split2.pfm/Results.ppp/Maps/Group.dscalar.nii',subs(s)));
    split1 = dscalar2double(split1,1); split2 = dscalar2double(split2,1);
    grp = ft_read_cifti(sprintf('PROFUMO/0All.pfm/Results.ppp/Maps/sub-%d.dscalar.nii',subs(s)));
    grp = dscalar2double(grp,1);
    assign = munkres(1-corr(split1,grp)); [i,~] = find(assign); split1 = split1(:,i); 
    A = corr(split1,grp); trt_grpcorrs(s,1:20) = A(eye(20)==1);
    assign = munkres(1-corr(split2,grp)); [i,~] = find(assign); split2 = split2(:,i); 
    A = corr(split2,grp); trt_grpcorrs(s+size(subs,2),1:20) = A(eye(20)==1);
    trt_grplabels(s,:) = cellstr(sprintf('sub-%d',subs(s)));
    
    A = corr(split1,split2); trt_corrs(s,:) = A(eye(20)==1);
    trt_labels(s,:) = cellstr(sprintf('sub-%d',subs(s)));
    clear split1 split2 grp A assign i
end

% Keep maps with median test-retest reliability > 0.4
corr_median = nanmedian(trt_corrs);
grp_median = nanmedian(trt_grpcorrs);
PFM_keep_group = find(corr_median>0.7 & grp_median>0.7);
length(PFM_keep_group)

example = ft_read_cifti('OLD/Overlap.dtseries.dtseries.nii'); example.dtseries = repmat(example.dtseries(:,1),1,length(PFM_keep_group));
notnans = find(isnan(example.dtseries(:,1))==0); arenans = find(isnan(example.dtseries(:,1))==1);
subject_maps_separate_runs = zeros(size(example.dtseries,1),length(PFM_keep_group)*size(subs,2));
subject_maps_group_run = zeros(size(example.dtseries,1),length(PFM_keep_group)*size(subs,2));
subject_maps_separate_split1 = zeros(size(example.dtseries,1),length(PFM_keep_group)*size(subs,2));
subject_maps_separate_split2 = zeros(size(example.dtseries,1),length(PFM_keep_group)*size(subs,2));
group_maps_group_run = zeros(size(example.dtseries,1),length(PFM_keep_group));
group_maps_group12_run = zeros(size(example.dtseries,1),length(PFM_keep_group));
PFM_keep_subjects = zeros(length(PFM_keep_group),size(subs,2)); PFM_keep_subjects_split1 = zeros(length(PFM_keep_group),size(subs,2)); PFM_keep_subjects_split2 = zeros(length(PFM_keep_group),size(subs,2));
corrs = zeros(2,length(PFM_keep_group)*size(subs,2));
missing = nan(length(PFM_keep_group),size(subs,2));
meanRs = zeros(length(subs),2);
for s = 1:size(subs,2)
    fprintf('Running subject %d (%d)\n',s,subs(s))
    subj = ft_read_cifti(sprintf('PROFUMO/2022_May/%d.pfm/Results.ppp/Maps/Group.dscalar.nii',subs(s)));
    grp = ft_read_cifti(sprintf('PROFUMO/0All.pfm/Results.ppp/Maps/sub-%d.dscalar.nii',subs(s)));
    split1 = ft_read_cifti(sprintf('PROFUMO/2022_May/%d_split1.pfm/Results.ppp/Maps/Group.dscalar.nii',subs(s)));
    split2 = ft_read_cifti(sprintf('PROFUMO/2022_May/%d_split2.pfm/Results.ppp/Maps/Group.dscalar.nii',subs(s)));
    subj = dscalar2double(subj,1); grp = dscalar2double(grp,1); split1 = dscalar2double(split1,1); split2 = dscalar2double(split2,1);
    assign = munkres(1-corr(subj,grp)); [i,~] = find(assign); subj = subj(:,i); PFM_keep_subjects(:,s) = i(PFM_keep_group); clear i assign
    assign = munkres(1-corr(split1,grp)); [i,~] = find(assign); split1 = split1(:,i); PFM_keep_subjects_split1(:,s) = i(PFM_keep_group); clear i assign
    assign = munkres(1-corr(split2,grp)); [i,~] = find(assign); split2 = split2(:,i); PFM_keep_subjects_split2(:,s) = i(PFM_keep_group); clear i assign
    grp = grp(:,PFM_keep_group); subj = subj(:,PFM_keep_group); split1 = split1(:,PFM_keep_group); split2 = split2(:,PFM_keep_group);
    
    R1 = corr(split1,split2); R1 = R1(eye(length(PFM_keep_group))==1); 
    R2 = corr(subj,grp); R2 = R2(eye(length(PFM_keep_group))==1);
    M = find(R1<0.2 | R2<0.2); missing(1:length(M),s) = M;
    grp(:,M) = nan; subj(:,M) = nan; split1(:,M) = nan; split2(:,M) = nan;
    R1(M) = nan; R2(M) = nan;
    corrs(1:2,(s-1)*length(PFM_keep_group)+1:s*length(PFM_keep_group)) = [R1 R2]';
    meanRs(s,:) = [nanmean(R1) nanmean(R2)];

    example.dtseries(notnans,:) = subj; example.hdr.dim(7) = length(PFM_keep_group); example.time = 1:length(PFM_keep_group);
    ft_write_cifti(sprintf('Results/Maps/Example_Subject_%d',subs(s)),example,'parameter','dtseries');

    subject_maps_separate_runs(notnans,(s-1)*length(PFM_keep_group)+1:s*length(PFM_keep_group)) = subj;
    subject_maps_group_run(notnans,(s-1)*length(PFM_keep_group)+1:s*length(PFM_keep_group)) = grp;
    subject_maps_separate_split1(notnans,(s-1)*length(PFM_keep_group)+1:s*length(PFM_keep_group)) = split1;
    subject_maps_separate_split2(notnans,(s-1)*length(PFM_keep_group)+1:s*length(PFM_keep_group)) = split2;
    
    clear subj grp split1 split2 R1 R2 M
end
grp = ft_read_cifti('PROFUMO/0All.pfm/Results.ppp/Maps/Group.dscalar.nii');
grp = dscalar2double(grp,1); grp = grp(:,PFM_keep_group);
grp12 = ft_read_cifti('PROFUMO/2022_May/mixedSubs.pfm/Results.ppp/Maps/Group.dscalar.nii');
grp12 = dscalar2double(grp12,1); 
assign = munkres(1-corr(grp12,grp)); [i,~] = find(assign); grp12 = grp12(:,i); clear i assign
group_maps_group_run(notnans,:) = grp;
group_maps_group12_run(notnans,:) = grp12;

% Run overlap script
[MinOverlap,ModeOverlap] = Overlap_tests(subject_maps_separate_runs,missing,subs,PFM_keep_group);

% Write out ciftis
example.dtseries = subject_maps_separate_runs; example.hdr.dim(7) = size(subject_maps_separate_runs,2); example.time = 1:size(subject_maps_separate_runs,2);
ft_write_cifti('Results/Maps/subject_maps_separate_runs',example,'parameter','dtseries');
example.dtseries = subject_maps_group_run; example.hdr.dim(7) = size(subject_maps_group_run,2); example.time = 1:size(subject_maps_group_run,2);
ft_write_cifti('Results/Maps/subject_maps_group_run',example,'parameter','dtseries');
example.dtseries = group_maps_group12_run; example.hdr.dim(7) = size(group_maps_group12_run,2); example.time = 1:size(group_maps_group12_run,2);
ft_write_cifti('Results/Maps/subject_maps_group12_run',example,'parameter','dtseries');
example.dtseries = subject_maps_separate_runs(:,(MinOverlap-1)*length(PFM_keep_group)+1:MinOverlap*length(PFM_keep_group)); example.hdr.dim(7) = length(PFM_keep_group); example.time = 1:length(PFM_keep_group);
ft_write_cifti(sprintf('Results/Maps/Example_Subject_Green_%d',subs(MinOverlap)),example,'parameter','dtseries');
example.dtseries = subject_maps_separate_runs(:,(ModeOverlap-1)*length(PFM_keep_group)+1:ModeOverlap*length(PFM_keep_group)); example.hdr.dim(7) = length(PFM_keep_group); example.time = 1:length(PFM_keep_group);
ft_write_cifti(sprintf('Results/Maps/Example_Subject_Red_%d',subs(ModeOverlap)),example,'parameter','dtseries');
example.dtseries = subject_maps_group_run(:,(MinOverlap-1)*length(PFM_keep_group)+1:MinOverlap*length(PFM_keep_group)); example.hdr.dim(7) = length(PFM_keep_group); example.time = 1:length(PFM_keep_group);
ft_write_cifti(sprintf('Results/Maps/Example_Subject_Green_%d_from_group',subs(MinOverlap)),example,'parameter','dtseries');
example.dtseries = subject_maps_group_run(:,(ModeOverlap-1)*length(PFM_keep_group)+1:ModeOverlap*length(PFM_keep_group)); example.hdr.dim(7) = length(PFM_keep_group); example.time = 1:length(PFM_keep_group);
ft_write_cifti(sprintf('Results/Maps/Example_Subject_Red_%d_from_group',subs(ModeOverlap)),example,'parameter','dtseries');
example.dtseries = subject_maps_separate_split1; example.hdr.dim(7) = size(subject_maps_separate_split1,2); example.time = 1:size(subject_maps_separate_split1,2);
ft_write_cifti('Results/Maps/subject_maps_separate_split1',example,'parameter','dtseries');
example.dtseries = subject_maps_separate_split2; example.hdr.dim(7) = size(subject_maps_separate_split2,2); example.time = 1:size(subject_maps_separate_split2,2);
ft_write_cifti('Results/Maps/subject_maps_separate_split2',example,'parameter','dtseries');
example.dtseries = group_maps_group_run; example.hdr.dim(7) = size(group_maps_group_run,2); example.time = 1:size(group_maps_group_run,2);
ft_write_cifti('Results/Maps/group_maps_group_run',example,'parameter','dtseries');
clear example 
save('PFM_keep.mat','PFM_keep_subjects','PFM_keep_group','PFM_keep_subjects_split1','PFM_keep_subjects_split2','subs','corrs','missing','meanRs','MinOverlap','ModeOverlap','arenans');

% Calculate spatial netmat
subject_maps_separate_runs(arenans,:) = []; subject_maps_group_run(arenans,:) = []; subject_maps_separate_split1(arenans,:) = []; subject_maps_separate_split2(arenans,:) = [];
Snet_subject_maps_separate_runs = zeros(length(PFM_keep_group),length(PFM_keep_group),size(subs,2));
Snet_subject_maps_group_run = zeros(length(PFM_keep_group),length(PFM_keep_group),size(subs,2));
Snet_subject_maps_separate_split1 = zeros(length(PFM_keep_group),length(PFM_keep_group),size(subs,2));
Snet_subject_maps_separate_split2 = zeros(length(PFM_keep_group),length(PFM_keep_group),size(subs,2));
for s = 1:size(subs,2)
    Snet_subject_maps_separate_runs(:,:,s) = corr(subject_maps_separate_runs(:,(s-1)*length(PFM_keep_group)+1:s*length(PFM_keep_group)),'rows','pairwise');
    Snet_subject_maps_group_run(:,:,s) = corr(subject_maps_group_run(:,(s-1)*length(PFM_keep_group)+1:s*length(PFM_keep_group)),'rows','pairwise');
    Snet_subject_maps_separate_split1(:,:,s) = corr(subject_maps_separate_split1(:,(s-1)*length(PFM_keep_group)+1:s*length(PFM_keep_group)),'rows','pairwise');
    Snet_subject_maps_separate_split2(:,:,s) = corr(subject_maps_separate_split2(:,(s-1)*length(PFM_keep_group)+1:s*length(PFM_keep_group)),'rows','pairwise');
end
save('Snets.mat','Snet_subject_maps_separate_runs','Snet_subject_maps_group_run','Snet_subject_maps_separate_split1','Snet_subject_maps_separate_split2')

% Calculate twin spatial correlations
Twins = [3 12; 4 16; 5 19; 7 18; 8 10; 11 17; 13 15; 14 20];
corrs_twins = zeros(2,8*length(PFM_keep_group)); in = 1;
for s = 1:size(Twins,1)
    s1a = subject_maps_separate_runs(:,(Twins(s,1)-1)*length(PFM_keep_group)+1 : Twins(s,1)*length(PFM_keep_group));
    s2a = subject_maps_separate_runs(:,(Twins(s,2)-1)*length(PFM_keep_group)+1 : Twins(s,2)*length(PFM_keep_group));
    s1b = subject_maps_group_run(:,(Twins(s,1)-1)*length(PFM_keep_group)+1 : Twins(s,1)*length(PFM_keep_group));
    s2b = subject_maps_group_run(:,(Twins(s,2)-1)*length(PFM_keep_group)+1 : Twins(s,2)*length(PFM_keep_group));
    for m = 1:length(PFM_keep_group)
        corrs_twins(1,in) = corr(s1a(:,m),s2a(:,m),'rows','pairwise');
        corrs_twins(2,in) = corr(s1b(:,m),s2b(:,m),'rows','pairwise'); in = in+1;
    end
end

% Calculate between-subject spatial correlations
NonTwins = 1:length(PFM_keep_group):size(subject_maps_group_run,2);
NonTwins(Twins(:,2)) = [];
I = ones(length(NonTwins)); I = triu(I,1); I = find(I==1);
corrs_between = zeros(2,length(I)*length(PFM_keep_group));
for m = 1:length(PFM_keep_group)
    A = corr(squeeze(subject_maps_separate_runs(:,NonTwins+m-1)),'rows','pairwise');
    corrs_between(1,(m-1)*length(I)+1:m*length(I)) = A(I); clear A
    A = corr(squeeze(subject_maps_group_run(:,NonTwins+m-1)),'rows','pairwise');
    corrs_between(2,(m-1)*length(I)+1:m*length(I)) = A(I); clear A
end

figure; set(gcf,'Position',[10 10 750 350],'PaperPositionMode','auto')
C = zeros(size(corrs,2),3); 
C((MinOverlap-1)*length(PFM_keep_group)+1:MinOverlap*length(PFM_keep_group),:) = repmat([0 1 0],length(PFM_keep_group),1); 
C((ModeOverlap-1)*length(PFM_keep_group)+1:ModeOverlap*length(PFM_keep_group),:) = repmat([1 0 0],length(PFM_keep_group),1); 
S = ones(size(corrs,2),1); S = S+10; 
S((MinOverlap-1)*length(PFM_keep_group)+1:MinOverlap*length(PFM_keep_group)) = 20;
S((ModeOverlap-1)*length(PFM_keep_group)+1:ModeOverlap*length(PFM_keep_group)) = 20;
swarmchart(ones(1,size(corrs,2)),corrs(1,:),S,C,'filled'); hold on; 
swarmchart(2*ones(1,size(corrs,2)),corrs(2,:),S,C,'filled'); 
swarmchart(3*ones(1,size(corrs_twins,2)),corrs_twins(1,:),10,'k','filled');
swarmchart(4*ones(1,size(corrs_twins,2)),corrs_twins(2,:),10,'k','filled');
swarmchart(5*ones(1,size(corrs_between,2)),corrs_between(1,:),10,'k','filled'); 
swarmchart(6*ones(1,size(corrs_between,2)),corrs_between(2,:),10,'k','filled'); 
set(gca,'xtick',[1 2 3 4 5 6],'xticklabel',{'test-retest reliability (6-scans vs 6-scans)','within-subject similarity (group-individual run)','twin similarity (individual runs)','twin similarity (group run)','between-subject similarity (individual runs)','between-subject similarity (group run)'},'FontSize',14)
ylabel('Spatial correlation','FontSize',14)
title('Spatial network similarity','FontSize',14)
axis([0.5 6.5 -0.4 1])
vline([1.5 2.5 3.5 4.5 5.5],':k')
print(gcf,'Results/swarm_similarity','-dpng','-r300');
MeanSTD = [nanmean(corrs(1,:)) nanstd(corrs(1,:)); ...
    nanmean(corrs(2,:)) nanstd(corrs(2,:));...
    nanmean(corrs_twins(1,:)) nanstd(corrs_twins(1,:));...
    nanmean(corrs_twins(2,:)) nanstd(corrs_twins(2,:));...
    nanmean(corrs_between(1,:)) nanstd(corrs_between(1,:));...
    nanmean(corrs_between(2,:)) nanstd(corrs_between(2,:))]
 
missing_matrix = zeros(size(missing));
for s = 1:length(subs)
    M = missing(:,s); M(isnan(M)) = []; 
    missing_matrix(M,s) = 1; clear M
end
figure; set(gcf,'Position',[615   301   739   246],'PaperPositionMode','auto')
imagesc(missing_matrix); colormap Bone
rectangle('Position',[MinOverlap-0.5 0.5 1 12],'EdgeColor','g')
rectangle('Position',[ModeOverlap-0.5 0.5 1 12],'EdgeColor','r')
ModeNames = {'L-FPN 1 (control 1)','ON 1 (visual 1)','D-FPN 1 (attention 1)','PN (somatosensory)','ON 2 (visual 2)','ON 3 (visual 4)','M-FPN 1 (default 1)','M-FPN 2 (default 2)','L-FPN 2 (control 2)','L-FPN 3 (control 3)','D-FPN 2 (attention 2)','L-FPN 4 (control 4)'};
set(gca,'xtick',1:1:length(subs),'xticklabel',subs,'ytick',1:12,'yticklabel',ModeNames,'FontSize',14)
ylabel('PFM modes','FontSize',14); xlabel('Subjecs','FontSize',14)
title('Number of missing modes (white = missing)','FontSize',14)
print(gcf,'Results/missing','-dpng','-r300');

% figure; 
% swarmchart(ones(1,size(subject_maps_group_run,1)*size(subject_maps_group_run,2)),subject_maps_group_run(:),10,'filled'); hold on; 
% swarmchart(2*ones(1,size(subject_maps_separate_runs,1)*size(subject_maps_separate_runs,2)),subject_maps_separate_runs(:),10,'filled'); 
% set(gca,'xtick',[1 2],'xticklabel',{'Subject maps from Group run','Subject maps from separate runs'});
% ylabel('PFM vertex values')
% title('PFM spatial weights')
% print(gcf,'Results/swarm_vertexvalues','-dpng','-r300');


clear all; close all; clc

load PFM_keep.mat
addpath('~/Dropbox/Matlab/hline_vline')

Tnet_subject_maps_separate_runs = zeros(length(PFM_keep_group),length(PFM_keep_group),size(subs,2));
Tnet_subject_maps_separate_split1 = zeros(length(PFM_keep_group),length(PFM_keep_group),size(subs,2));
Tnet_subject_maps_separate_split2 = zeros(length(PFM_keep_group),length(PFM_keep_group),size(subs,2));
Tnet_subject_maps_group_run = zeros(length(PFM_keep_group),length(PFM_keep_group),size(subs,2));

for s = 1:size(subs,2)
    fprintf('Running subject %d (%d)\n',s,subs(s))
    subj = load(sprintf('PROFUMO/2022_May/%d.pfm/Results.ppp/NetMats/Group.csv',subs(s)));
    split1 = load(sprintf('PROFUMO/2022_May/%d_split1.pfm/Results.ppp/NetMats/Group.csv',subs(s)));
    split2 = load(sprintf('PROFUMO/2022_May/%d_split2.pfm/Results.ppp/NetMats/Group.csv',subs(s)));
    grp = load(sprintf('PROFUMO/0All.pfm/Results.ppp/Netmats/sub-%d.csv',subs(s)));
    M = missing(:,s); M(isnan(M)) = [];
    A = grp(PFM_keep_group,PFM_keep_group); A(M,:) = nan; A(:,M) = nan;
    Tnet_subject_maps_group_run(:,:,s) = A; clear A
    A = subj(PFM_keep_subjects(:,s),PFM_keep_subjects(:,s)); A(M,:) = nan; A(:,M) = nan;
    Tnet_subject_maps_separate_runs(:,:,s) = A; clear A
    A = split1(PFM_keep_subjects_split1(:,s),PFM_keep_subjects_split1(:,s));  A(M,:) = nan; A(:,M) = nan;
    Tnet_subject_maps_separate_split1(:,:,s) = A; clear A
    A = split2(PFM_keep_subjects_split2(:,s),PFM_keep_subjects_split2(:,s)); A(M,:) = nan; A(:,M) = nan;
    Tnet_subject_maps_separate_split2(:,:,s) = A; clear A
clear subj split1 split2 grp
end
grp = load('PROFUMO/0All.pfm/Results.ppp/Netmats/Group.csv'); grp = grp(PFM_keep_group,PFM_keep_group);
Tnet_subject_maps_group_run(Tnet_subject_maps_group_run==-1) = NaN;
Tnet_subject_maps_separate_runs(Tnet_subject_maps_separate_runs==-1) = NaN;
grp(grp==-1) = NaN;

% R to z transform
load Snets.mat
Tnet_subject_maps_group_run = 0.5*log((1+Tnet_subject_maps_group_run)./(1-Tnet_subject_maps_group_run));
Tnet_subject_maps_separate_runs = 0.5*log((1+Tnet_subject_maps_separate_runs)./(1-Tnet_subject_maps_separate_runs));
Tnet_subject_maps_separate_split1 = 0.5*log((1+Tnet_subject_maps_separate_split1)./(1-Tnet_subject_maps_separate_split1));
Tnet_subject_maps_separate_split1 = 0.5*log((1+Tnet_subject_maps_separate_split1)./(1-Tnet_subject_maps_separate_split1));
Snet_subject_maps_group_run = 0.5*log((1+Snet_subject_maps_group_run)./(1-Snet_subject_maps_group_run));
Snet_subject_maps_separate_runs = 0.5*log((1+Snet_subject_maps_separate_runs)./(1-Snet_subject_maps_separate_runs));
Snet_subject_maps_separate_split1 = 0.5*log((1+Snet_subject_maps_separate_split1)./(1-Snet_subject_maps_separate_split1));
Snet_subject_maps_separate_split1 = 0.5*log((1+Snet_subject_maps_separate_split1)./(1-Snet_subject_maps_separate_split1));

% Figures
Snet_subject_maps_group_run(Snet_subject_maps_group_run==1) = NaN;
Snet_subject_maps_separate_runs(Snet_subject_maps_separate_runs==1) = NaN;
for s = [MinOverlap ModeOverlap]
    if s==MinOverlap; C = 'g'; elseif s==ModeOverlap; C = 'r'; end
    figure; set(gcf,'Position',[10 10 900 800],'PaperPositionMode','auto')
    minmax = [-0.45 0.45];
    subplot(3,2,1); imagesc(nanmean(Snet_subject_maps_group_run,3),minmax); colorbar; title('Snet Group','FontSize',14)
    subplot(3,2,2); imagesc(grp,minmax); colorbar; title('Tnet group','FontSize',12)
    subplot(3,2,3); imagesc(Snet_subject_maps_group_run(:,:,s),minmax); colorbar; title(sprintf('Snet subject %d from group run',subs(s)),'Color',C,'FontSize',14)
    subplot(3,2,4); imagesc(Tnet_subject_maps_group_run(:,:,s),minmax); colorbar; title(sprintf('Tnet subject %d from group run',subs(s)),'Color',C,'FontSize',14)
    subplot(3,2,5); imagesc(Snet_subject_maps_separate_runs(:,:,s),minmax); colorbar; title(sprintf('Snet subject %d from separate run',subs(s)),'Color',C,'FontSize',14)
    subplot(3,2,6); imagesc(Tnet_subject_maps_separate_runs(:,:,s),minmax); colorbar; title(sprintf('Tnet subject %d from separate run',subs(s)),'Color',C,'FontSize',14)
    if s==MinOverlap
        print(gcf,sprintf('Results/TnetSnet_Example_Subject_Green_%d',subs(s)),'-dpng','-r300');
    end
    if s==ModeOverlap
        print(gcf,sprintf('Results/TnetSnet_Example_Subject_Red_%d',subs(s)),'-dpng','-r300');
    end
end

% Test - retest figure for tnet and snet
I = ones(length(PFM_keep_group)); I = triu(I,1); I = find(I==1);
Snet_grp = zeros(length(I),length(subs));
Snet_split1 = zeros(length(I),length(subs));
Snet_split2 = zeros(length(I),length(subs));
Snet_subj = zeros(length(I),length(subs));
Tnet_grp = zeros(length(I),length(subs));
Tnet_split1 = zeros(length(I),length(subs));
Tnet_split2 = zeros(length(I),length(subs));
Tnet_subj = zeros(length(I),length(subs));
Snet_corrs = zeros(2,size(subs,2));
Tnet_corrs = zeros(2,size(subs,2));
for s = 1:length(subs)
    A = squeeze(Snet_subject_maps_separate_runs(:,:,s));
    Snet_subj(:,s) = A(I); clear A
    A = squeeze(Snet_subject_maps_separate_split1(:,:,s));
    Snet_split1(:,s) = A(I); clear A
    A = squeeze(Snet_subject_maps_separate_split2(:,:,s));
    Snet_split2(:,s) = A(I); clear A
    A = squeeze(Snet_subject_maps_group_run(:,:,s));
    Snet_grp(:,s) = A(I); clear A
    A = squeeze(Tnet_subject_maps_separate_split1(:,:,s));
    Tnet_split1(:,s) = A(I); clear A
    A = squeeze(Tnet_subject_maps_separate_split2(:,:,s));
    Tnet_split2(:,s) = A(I); clear A
    A = squeeze(Tnet_subject_maps_separate_runs(:,:,s));
    Tnet_subj(:,s) = A(I); clear A
    A = squeeze(Tnet_subject_maps_group_run(:,:,s));
    Tnet_grp(:,s) = A(I); clear A

    R1 = corr(Tnet_split1(:,s),Tnet_split2(:,s),'rows','pairwise'); 
    R2 = corr(Tnet_subj(:,s),Tnet_grp(:,s),'rows','pairwise'); 
    Tnet_corrs(1:2,s) = [R1; R2];

    R1 = corr(Snet_split1(:,s),Snet_split2(:,s),'rows','pairwise'); 
    R2 = corr(Snet_subj(:,s),Snet_grp(:,s),'rows','pairwise');
    Snet_corrs(1:2,s) = [R1; R2];     
end

% Calculate twin Tnet and Snet correlations
Twins = [3 12; 4 16; 5 19; 7 18; 8 10; 11 17; 13 15; 14 20];
Tnet_corrs_twins = zeros(2,8); Snet_corrs_twins = zeros(2,8); in = 1;
for s = 1:size(Twins,1)
    Tnets1a = Tnet_subj(:,Twins(s,1));
    Tnets2a = Tnet_subj(:,Twins(s,2));
    Tnets1b = Tnet_grp(:,Twins(s,1));
    Tnets2b = Tnet_grp(:,Twins(s,2));
    Snets1a = Snet_subj(:,Twins(s,1));
    Snets2a = Snet_subj(:,Twins(s,2));
    Snets1b = Snet_grp(:,Twins(s,1));
    Snets2b = Snet_grp(:,Twins(s,2));
    Tnet_corrs_twins(1,s) = corr(Tnets1a,Tnets2a,'rows','pairwise');
    Tnet_corrs_twins(2,s) = corr(Tnets1b,Tnets2b,'rows','pairwise');
    Snet_corrs_twins(1,s) = corr(Snets1a,Snets2a,'rows','pairwise');
    Snet_corrs_twins(2,s) = corr(Snets1b,Snets2b,'rows','pairwise'); in = in+1;
end

% Calculate between-subject Tnet and Snet correlations
I = ones(12); I = triu(I,1); I = find(I==1);
Tnet_corrs_between = zeros(2,length(I)); Snet_corrs_between = zeros(2,length(I));
A = corr(Tnet_subj(:,setdiff(1:20,Twins(:,2))),'rows','pairwise'); Tnet_corrs_between(1,:) = A(I); clear A
A = corr(Snet_subj(:,setdiff(1:20,Twins(:,2))),'rows','pairwise'); Snet_corrs_between(1,:) = A(I); clear A
A = corr(Tnet_grp(:,setdiff(1:20,Twins(:,2))),'rows','pairwise'); Tnet_corrs_between(2,:) = A(I); clear A
A = corr(Snet_grp(:,setdiff(1:20,Twins(:,2))),'rows','pairwise'); Snet_corrs_between(2,:) = A(I); clear A

figure; set(gcf,'Position',[10 10 750 350],'PaperPositionMode','auto')
C = zeros(length(subs),3); C(MinOverlap,:) = [0 1 0]; C(ModeOverlap,:) = [1 0 0];
S = ones(length(subs),1); S = S+10; S(MinOverlap) = 20; S(ModeOverlap) = 20; 
swarmchart(ones(1,size(Tnet_corrs,2)),Tnet_corrs(1,:),S,C,'filled'); hold on; 
swarmchart(2*ones(1,size(Tnet_corrs,2)),Tnet_corrs(2,:),S,C,'filled'); 
swarmchart(3*ones(1,size(Tnet_corrs_twins,2)),Tnet_corrs_twins(1,:),10,'k','filled'); 
swarmchart(4*ones(1,size(Tnet_corrs_twins,2)),Tnet_corrs_twins(2,:),10,'k','filled'); 
swarmchart(5*ones(1,size(Tnet_corrs_between,2)),Tnet_corrs_between(1,:),10,'k','filled'); 
swarmchart(6*ones(1,size(Tnet_corrs_between,2)),Tnet_corrs_between(2,:),10,'k','filled'); 
set(gca,'xtick',[1 2 3 4 5 6],'xticklabel',{'test-retest reliability','within-subject similarity (group-individual run)','twin similarity (individual runs)','twin similarity (group run)','between-subject similarity (individual runs)','between-subject similarity (group run)'},'FontSize',14)
ylabel('Correlation','FontSize',14)
title('Temporal netmat similarity','FontSize',14)
axis([0.5 6.5 -0.4 1])
vline([1.5 2.5 3.5 4.5 5.5],':k')
print(gcf,'Results/swarm_similarity_Tnets','-dpng','-r300');


figure; set(gcf,'Position',[10 10 750 350],'PaperPositionMode','auto')
swarmchart(ones(1,size(Snet_corrs,2)),Snet_corrs(1,:),S,C,'filled'); hold on; 
swarmchart(2*ones(1,size(Snet_corrs,2)),Snet_corrs(2,:),S,C,'filled'); 
swarmchart(3*ones(1,size(Snet_corrs_twins,2)),Snet_corrs_twins(1,:),10,'k','filled'); 
swarmchart(4*ones(1,size(Snet_corrs_twins,2)),Snet_corrs_twins(2,:),10,'k','filled');
swarmchart(5*ones(1,size(Snet_corrs_between,2)),Snet_corrs_between(1,:),10,'k','filled'); 
swarmchart(6*ones(1,size(Snet_corrs_between,2)),Snet_corrs_between(2,:),10,'k','filled'); 
set(gca,'xtick',[1 2 3 4 5 6],'xticklabel',{'test-retest reliability','within-subject similarity (group-individual run)','twin similarity (individual runs)','twin similarity (group run)','between-subject similarity (individual runs)','between-subject similarity (group run)'},'FontSize',14)
ylabel('Correlation','FontSize',14)
title('Spatial netmat similarity','FontSize',14)
axis([0.5 6.5 -0.4 1])
vline([1.5 2.5 3.5 4.5 5.5],':k')
print(gcf,'Results/swarm_similarity_Snets','-dpng','-r300');




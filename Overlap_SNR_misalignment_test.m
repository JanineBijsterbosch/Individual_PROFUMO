clear all; close all; clc

D = dir('PROFUMO/2022_May/mixedSubs.pfm/Results.ppp/Maps/sub*');
load("PFM_keep.mat")

Overlap_volume = zeros(length(D),10,3);
I = ones(length(PFM_keep_group)); I = triu(I,1); I = find(I==1);
for s = 1:length(D)
    fprintf('Running subject %d\n',s)
    sub = D(s).name(5:10);
    mis = missing(:,subs==str2num(sub)); mis(isnan(mis)==1) = [];

    subj = ft_read_cifti(sprintf('PROFUMO/2022_May/%s.pfm/Results.ppp/Maps/Group.dscalar.nii',sub)); 
    subj = dscalar2double(subj,1); subj = subj(:,PFM_keep_subjects(:,subs==str2num(sub))); subj(:,mis) = [];
    M = subj; M(M<1) = 0; M(abs(M)>=1) = 1;
    M = nansum(M,2); M(M==1) = 0;
    for v = 1:size(Overlap_volume,2)
        Overlap_volume(s,v,1) = length(find(M==v));
    end
    
    grp = ft_read_cifti(sprintf('PROFUMO/0all.pfm/Results.ppp/Maps/sub-%s.dscalar.nii',sub)); 
    grp = dscalar2double(grp,1); grp = grp(:,PFM_keep_group); grp(:,mis) = [];
    grp = grp./repmat(max(grp),size(grp,1),1); grp = grp.*repmat(max(subj),size(grp,1),1); 
    M = grp; M(M<1) = 0; M(abs(M)>=1) = 1;
    M = nansum(M,2); M(M==1) = 0;
    for v = 1:size(Overlap_volume,2)
        Overlap_volume(s,v,2) = length(find(M==v));
    end

    mixed = ft_read_cifti(sprintf('PROFUMO/2022_May/mixedSubs.pfm/Results.ppp/Maps/%s',D(s).name)); 
    mixed = dscalar2double(mixed,1);
    assign = munkres(1-corr(mixed,grp)); [i,~] = find(assign); mixed = mixed(:,i);
    mixed = mixed./repmat(max(mixed),size(mixed,1),1); mixed = mixed.*repmat(max(subj),size(mixed,1),1); 
    M = mixed; M(M<1) = 0; M(abs(M)>=1) = 1;
    M = nansum(M,2); M(M==1) = 0;
    for v = 1:size(Overlap_volume,2)
        Overlap_volume(s,v,3) = length(find(M==v));
    end

end

figure; set(gcf,'Position',[600 600 700 200],'PaperPositionMode','auto')
swarmchart(ones(1,size(Overlap_volume,1)),squeeze(Overlap_volume(:,2,1)),10,[0.9290 0.6940 0.1250],'filled'); hold on; 
swarmchart(2*ones(1,size(Overlap_volume,1)),squeeze(Overlap_volume(:,2,2)),10,[0.4940 0.1840 0.5560],'filled'); hold on; 
swarmchart(3*ones(1,size(Overlap_volume,1)),squeeze(Overlap_volume(:,2,3)),10,[0.8500 0.3250 0.0980],'filled'); hold on; 

swarmchart(5*ones(1,size(Overlap_volume,1)),squeeze(Overlap_volume(:,3,1)),10,[0.9290 0.6940 0.1250],'filled'); hold on; 
swarmchart(6*ones(1,size(Overlap_volume,1)),squeeze(Overlap_volume(:,3,2)),10,[0.4940 0.1840 0.5560],'filled'); hold on; 
swarmchart(7*ones(1,size(Overlap_volume,1)),squeeze(Overlap_volume(:,3,3)),10,[0.8500 0.3250 0.0980],'filled'); hold on; 

swarmchart(9*ones(1,size(Overlap_volume,1)),squeeze(Overlap_volume(:,4,1)),10,[0.9290 0.6940 0.1250],'filled'); hold on; 
swarmchart(10*ones(1,size(Overlap_volume,1)),squeeze(Overlap_volume(:,4,2)),10,[0.4940 0.1840 0.5560],'filled'); hold on; 
swarmchart(11*ones(1,size(Overlap_volume,1)),squeeze(Overlap_volume(:,4,3)),10,[0.8500 0.3250 0.0980],'filled'); hold on; 

set(gca,'xtick',[2 6 10],'xticklabel',{'two maps overlap','three maps overlap','four maps overlap'})
legend({'Individual subject PROFUMO (12 runs from 1 participant)','Group PROFUMO (12 runs from each of 20 participants)','Group PROFUMO matched data size (1 run from each of 12 participants)'})
ylabel('Number of vertices')
title('Spatial network overlap comparison between PROFUMO runs')
print(gcf,'Results/swarm_overlap_SNR_misalignment_test','-dpng','-r300');





function [MinOverlap,ModeOverlap] = Overlap_tests(subj_maps, missing, subs, PFM_keep_group)

addpath(genpath('~/Dropbox/Matlab/fieldtrip-20190819/'),'-END')
example = ft_read_cifti('OLD/Overlap.dtseries.dtseries.nii');
notnans = find(isnan(example.dtseries(:,1))==0); arenans = find(isnan(example.dtseries(:,1))==1);
subj_maps(arenans,:) = [];

Overlap_volume = zeros(length(subs),10);
FTB = zeros(size(example.dtseries,1),length(subs));
Overlap_maps = zeros(size(example.dtseries,1),length(subs));
spatial_overlap_matrix = zeros(length(PFM_keep_group),length(PFM_keep_group),length(subs));
for s = 1:length(subs)
    M = subj_maps(:,(s-1)*length(PFM_keep_group)+1:s*length(PFM_keep_group));
    [w,I] = max(M,[],2); I(w<1) = nan;
    FTB(notnans,s) = I;
    spatial_overlap_matrix(:,:,s) = corr(M);
    spatial_overlap_matrix(:,:,s) = 0.5*log((1+spatial_overlap_matrix(:,:,s))./(1-spatial_overlap_matrix(:,:,s)));
    M(M<1) = 0; M(abs(M)>=1) = 1;
    M = nansum(M,2); M(M==1) = 0;
    Overlap_maps(notnans,s) = M;
    for v = 1:size(Overlap_volume,2)
        Overlap_volume(s,v) = length(find(M==v));
    end
end
clear v s M 

example.dtseries = Overlap_maps; example.hdr.dim(7) = size(Overlap_maps,2); example.time = 1:size(Overlap_maps,2);
ft_write_cifti('Results/Maps/overlap_maps',example,'parameter','dtseries');
example.dtseries = FTB; 
ft_write_cifti('Results/Maps/Find_the_biggest',example,'parameter','dtseries');
clear example notnans arenans

% Select subjects
NoMissing = find(nansum(missing,1)==0);
[~,MinOverlap] = min(Overlap_volume(NoMissing,2)); MinOverlap = NoMissing(MinOverlap);
[~,ModeOverlap] = mode(Overlap_volume(NoMissing,2)); ModeOverlap = NoMissing(ModeOverlap);

% Percentile of threshold (1)
T = length(find(subj_maps(:)>=1))/(91282*240)*100;
fprintf('Spatial threshold of 1 is %1.1f percentile\n',T);

figure; set(gcf,'Position',[10 10 900 700],'PaperPositionMode','auto')
subplot(2,2,1);
imagesc(squeeze(nanmean(spatial_overlap_matrix,3)),[-0.5 0.5]);
ModeNames = {'L-FPN 1 (control 1)','ON 1 (visual 1)','D-FPN 1 (attention 1)','PN (somatosensory)','ON 2 (visual 2)','ON 3 (visual 4)','M-FPN 1 (default 1)','M-FPN 2 (default 2)','L-FPN 2 (control 2)','L-FPN 3 (control 3)','D-FPN 2 (attention 2)','L-FPN 4 (control 4)'};
set(gca,'xtick',1:length(PFM_keep_group),'xticklabel',ModeNames,'ytick',1:length(PFM_keep_group),'yticklabel',ModeNames,'FontSize',14)
a=colorbar; ylabel(a,'Mean Z-transformed correlation','Rotation',270,'FontSize',14); a.Label.Position(1) = 3;
title('A. Spatial overlap matrix averaged across participants','FontSize',14)

subplot(2,2,2);
C = zeros(length(subs),3); C(MinOverlap,:) = [0 1 0]; C(ModeOverlap,:) = [1 0 0];
S = ones(length(subs),1); S = S+10; S(MinOverlap) = 20; S(ModeOverlap) = 20; 
swarmchart(ones(1,size(Overlap_volume,1)),Overlap_volume(:,2),S,C,'filled'); hold on; 
swarmchart(2*ones(1,size(Overlap_volume,1)),Overlap_volume(:,3),S,C,'filled'); 
swarmchart(3*ones(1,size(Overlap_volume,1)),Overlap_volume(:,4),S,C,'filled');  
swarmchart(4*ones(1,size(Overlap_volume,1)),sum(Overlap_volume(:,5:end),2),S,C,'filled'); 
set(gca,'xtick',[1 2 3 4],'xticklabel',{'two maps overlap','three maps overlap','four maps overlap','five or more maps overlap'},'FontSize',14)
ylabel('Number of vertices','FontSize',14)
title('B. Spatial network overlap count','FontSize',14)

subplot(2,2,3);
A = imread(sprintf('Results/wb_view_screenshots/Overlap_sub%02d_new.png',MinOverlap));
imshow(A,'border','tight'); 
title(strcat('\color{black}C. Overlap map for participant \color{green}',sprintf('%d',subs(MinOverlap))),'FontSize',14);

subplot(2,2,4);
A = imread(sprintf('Results/wb_view_screenshots/Overlap_sub%02d_new.png',ModeOverlap));
imshow(A,'border','tight'); 
title(strcat('\color{black}D. Overlap map for participant \color{red}',sprintf('%d',subs(ModeOverlap))),'FontSize',14);

print(gcf,'Results/swarm_overlap','-dpng','-r300');



clear all; close all; clc

%% Run on CPHC
% addpath(genpath('/scratch/janine.bijsterbosch/HCP/Matlab/'),'-END')
% 
% subs = [105923 114823 125525 130518 137128 144226 146129 158035 169343 177746 185442 192439 195041 200614 204521 562345 601127 783462 859671 861456];
% 
% % Test-retest reliability and group-subject similarity
% trt_repeats = nan(size(subs,2)*2,20); 
% for s = 1:size(subs,2)
%     fprintf('Running subject %d (%d)\n',s,subs(s))
%     split1 = ft_read_cifti(sprintf('PROFUMO/2022_May/%d.pfm/Results.ppp/Maps/Group.dscalar.nii',subs(s)));
%     split2 = ft_read_cifti(sprintf('PROFUMO/2022_May/%d_repeat.pfm/Results.ppp/Maps/Group.dscalar.nii',subs(s)));
%     split1 = dscalar2double(split1,1); split2 = dscalar2double(split2,1);
%     grp = ft_read_cifti(sprintf('PROFUMO/0All.pfm/Results.ppp/Maps/sub-%d.dscalar.nii',subs(s)));
%     grp = dscalar2double(grp,1);
%     assign = munkres(1-corr(split1,grp)); [i,~] = find(assign); split1 = split1(:,i); 
%     assign = munkres(1-corr(split2,grp)); [i,~] = find(assign); split2 = split2(:,i);
%     A = corr(split1,split2); trt_repeats(s,1:20) = A(eye(20)==1);
% end
% 
% save('Subject_repeats.mat','trt_repeats');

%% Run locally after grapping Subject_repeats.mat from CHPC
load Subject_repeats.mat
load PFM_keep.mat PFM_keep_group missing
repeats12 = zeros(20,12);
for s = 1:20
    repeats12(s,:) = trt_repeats(s,PFM_keep_group);
    M = missing(:,s); M(isnan(M)==1) = [];
    repeats12(s,M) = nan;
end
repeats12 = reshape(repeats12,20*12,1);
figure
swarmchart(ones(1,size(repeats12,1)),repeats12,10,"blue","filled"); hold on
plot([0.75 1.25],[nanmean(repeats12) nanmean(repeats12)],'k'); 
plot([0.75 1.25],[nanmedian(repeats12) nanmedian(repeats12)],'k','linewidth',2);
title('Repeat PROFUMO similarity of subject maps')
print(gcf,'Results/October_feedback_Repeat_PROFUMO_similarity','-dpng','-r300');


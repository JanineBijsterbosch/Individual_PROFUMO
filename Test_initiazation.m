clear all; close all; clc

addpath(genpath('~/Dropbox/Matlab/fieldtrip-20190819/'),'-END')
addpath('~/Dropbox/Matlab/','-END')
addpath(genpath('~/Dropbox/Matlab/FSLNets'),'-END')
%addpath(genpath('~/Dropbox/Matlab/Violinplot-Matlab-master'))
addpath(genpath('~/Dropbox/Matlab/ICC'),'-END')

load PFM_keep.mat

% Relate group maps to initiation stages for group run:
IT_Rs = zeros(151,length(PFM_keep_group));
grp = ft_read_cifti('PROFUMO/0All.pfm/Results.ppp/Maps/Group.dscalar.nii');
grp = dscalar2double(grp,1); grp = grp(:,PFM_keep_group); IT = 1;
D = dir('PROFUMO/0All.pfm/Intermediates/Model1/GroupMaps*'); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/0All.pfm/Intermediates/Model1/%s',D(n).name),'/dataset');
    data = data(:,PFM_keep_group);
    R = corr(grp,data);
    IT_Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end
D = dir('PROFUMO/0All.pfm/Intermediates/Model2/GroupMaps*'); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/0All.pfm/Intermediates/Model2/%s',D(n).name),'/dataset');
    data = data(:,PFM_keep_group);
    R = corr(grp,data);
    IT_Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end
D = dir('PROFUMO/0All.pfm/Intermediates/Model3/GroupMaps*'); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/0All.pfm/Intermediates/Model3/%s',D(n).name),'/dataset');
    data = data(:,PFM_keep_group);
    R = corr(grp,data);
    IT_Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end
D = dir('PROFUMO/0All.pfm/Intermediates/Model4/GroupMaps*'); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/0All.pfm/Intermediates/Model4/%s',D(n).name),'/dataset');
    data = data(:,PFM_keep_group);
    R = corr(grp,data);
    IT_Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end
D = dir('PROFUMO/0All.pfm/Intermediates/Model5/GroupMaps*'); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/0All.pfm/Intermediates/Model5/%s',D(n).name),'/dataset');
    data = data(:,PFM_keep_group);
    R = corr(grp,data);
    IT_Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end

% Relate group maps to initiation stages for subject 'run':
IT_s1Rs = zeros(151,length(PFM_keep_group));
sub = subs(MinOverlap);
grp = ft_read_cifti(sprintf('PROFUMO/2022_May/%d.pfm/Results.ppp/Maps/Group.dscalar.nii',sub));
grp = dscalar2double(grp,1); grp = grp(:,PFM_keep_subjects(:,MinOverlap)); IT = 1;
D = dir(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model1/GroupMaps*',sub)); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model1/%s',sub,D(n).name),'/dataset');
    data = data(:,PFM_keep_subjects(:,MinOverlap));
    R = corr(grp,data);
    IT_s1Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end
D = dir(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model2/GroupMaps*',sub)); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model2/%s',sub,D(n).name),'/dataset');
    data = data(:,PFM_keep_subjects(:,MinOverlap));
    R = corr(grp,data);
    IT_s1Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end
D = dir(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model3/GroupMaps*',sub)); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model3/%s',sub,D(n).name),'/dataset');
    data = data(:,PFM_keep_subjects(:,MinOverlap));
    R = corr(grp,data);
    IT_s1Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end
D = dir(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model4/GroupMaps*',sub)); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model4/%s',sub,D(n).name),'/dataset');
    data = data(:,PFM_keep_subjects(:,MinOverlap));
    R = corr(grp,data);
    IT_s1Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end
D = dir(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model5/GroupMaps*',sub)); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model5/%s',sub,D(n).name),'/dataset');
    data = data(:,PFM_keep_subjects(:,MinOverlap));
    R = corr(grp,data);
    IT_s1Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end

% Relate group maps to initiation stages for subject 'run':
IT_s2Rs = zeros(151,length(PFM_keep_group));
sub = subs(ModeOverlap);
grp = ft_read_cifti(sprintf('PROFUMO/2022_May/%d.pfm/Results.ppp/Maps/Group.dscalar.nii',sub));
grp = dscalar2double(grp,1); grp = grp(:,PFM_keep_subjects(:,ModeOverlap)); IT = 1;
D = dir(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model1/GroupMaps*',sub)); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model1/%s',sub,D(n).name),'/dataset');
    data = data(:,PFM_keep_subjects(:,ModeOverlap));
    R = corr(grp,data);
    IT_s2Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end
D = dir(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model2/GroupMaps*',sub)); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model2/%s',sub,D(n).name),'/dataset');
    data = data(:,PFM_keep_subjects(:,ModeOverlap));
    R = corr(grp,data);
    IT_s2Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end
D = dir(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model3/GroupMaps*',sub)); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model3/%s',sub,D(n).name),'/dataset');
    data = data(:,PFM_keep_subjects(:,ModeOverlap));
    R = corr(grp,data);
    IT_s2Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end
D = dir(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model4/GroupMaps*',sub)); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model4/%s',sub,D(n).name),'/dataset');
    data = data(:,PFM_keep_subjects(:,ModeOverlap));
    R = corr(grp,data);
    IT_s2Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end
D = dir(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model5/GroupMaps*',sub)); 
for n = 1:length(D)
    data = h5read(sprintf('PROFUMO/2022_May/%d.pfm/Intermediates/Model5/%s',sub,D(n).name),'/dataset');
    data = data(:,PFM_keep_subjects(:,ModeOverlap));
    R = corr(grp,data);
    IT_s2Rs(IT,:) = R(eye(length(PFM_keep_group))==1); IT = IT+1;
end

figure; 
subplot(3,1,1); errorbar(mean(IT_Rs,2),std(IT_Rs,[],2)); axis([1 151 0.9 1])
title('Spatial correlation between final group maps and iterations for Traditional group PROFUMO')
subplot(3,1,2); errorbar(mean(IT_s1Rs,2),std(IT_s1Rs,[],2)); axis([1 151 0.9 1])
title(sprintf('Spatial correlation between final group maps and iterations for Example Green subject (sub %d)',subs(MinOverlap)),'Color','g')
subplot(3,1,3); errorbar(mean(IT_s2Rs,2),std(IT_s2Rs,[],2)); axis([1 151 0.9 1])
title(sprintf('Spatial correlation between final group maps and iterations for Example Red subject (sub %d)',subs(ModeOverlap)),'Color','r')
set(gcf,'Position',[620   457   731   590],'PaperPositionMode','auto')
print(gcf,'Results/iteration','-dpng','-r300');
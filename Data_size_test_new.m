clear all; close all; clc

What_to_do = 'plot'; % 'run' or 'plot'

if strcmp(What_to_do,'run')
    addpath(genpath('~/Dropbox/Matlab/fieldtrip-20190819/'),'-END')

    subs = [105923 114823 125525 130518 137128 144226 146129 158035 169343 177746 185442 192439 195041 200614 204521 562345 601127 783462 859671 861456];
    Gold_standard_maps = 'subject_run'; % 'group_run' or 'subject_run'
    data_size = [100 200 300 400 500 600 700 800 900 1000 1100];
    load PFM_keep.mat

    % Test-retest reliability and group-subject similarity
    corrs = nan(size(data_size,2),size(subs,2),12);
    for s = 1:size(subs,2)
        subj = ft_read_cifti(sprintf('Results/Maps/Example_Subject_%d.dtseries.nii',subs(s)));
        subj = subj.dtseries; subj(arenans,:) = []; subj(:,isnan(subj(1,:))==1) = [];
        %     grp = ft_read_cifti(sprintf('Results/Maps/Example_Subject_%d_from_group.dtseries.nii',subs(s)));
        %     grp = grp.dtseries; grp(grp(:,1)==0,:) = [];
        for d = 1:size(data_size,2)
            fprintf('Running subject %d (%d) with %d runs\n',s,subs(s),data_size(d));
            test = ft_read_cifti(sprintf('PROFUMO/2022_May/Data_amount_test_new/%d_%04d.pfm/Results.ppp/Maps/Group.dscalar.nii',subs(s),data_size(d)));
            test = dscalar2double(test,1);
            assign = munkres(1-corr(test,subj)); [i,~] = find(assign); test = test(:,i);
            if strcmp(Gold_standard_maps,'group_run')
                A = corr(test,grp); corrs(d,s,1:size(subj,2)) = A(eye(size(subj,2))==1);
            elseif strcmp(Gold_standard_maps,'subject_run')
                A = corr(test,subj); corrs(d,s,1:size(subj,2)) = A(eye(size(subj,2))==1);
            end
            clear split1 split2 grp A assign i
        end
    end

    corrs = reshape(corrs,size(data_size,2),size(subs,2)*12);
    save('Data_size_test_new.mat','corrs');

elseif strcmp(What_to_do,'plot')

    load Data_size_test_new.mat
    data_size = 100:100:1100;
    data_size_7T = 75:75:825;

    figure; set(gcf,'Position',[286 716 750 331],'PaperPositionMode','auto');
    for d = 1:length(data_size)
        swarmchart(ones(1,size(corrs,2))*d,corrs(d,:),10,'filled'); hold on;
        plot([0.75+(d-1) 1.25+(d-1)],[nanmean(corrs(d,:)) nanmean(corrs(d,:))],'k'); plot([0.75+(d-1) 1.25+(d-1)],[nanmedian(corrs(d,:)) nanmedian(corrs(d,:))],'k','linewidth',2);
    end
    set(gca,'xtick',1:12,'xticklabel',8*data_size + 4*data_size_7T,'FontSize',14)
    xlabel('Number of TRs included (split across all 12 runs)','FontSize',14)
    ylabel('Correlation with subject maps from all 12 runs','FontSize',14)
    axis([0.5 11.5 -0.2 1])
    title('Amount of data needed for individual PROFUMO','FontSize',14)
    print(gcf,'Results/data_amount_new','-dpng','-r300');
end

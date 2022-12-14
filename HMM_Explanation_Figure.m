clear all; close all; clc

addpath(genpath('~/Dropbox/Matlab/fieldtrip-20190819/'),'-END')

subs = 125525;

% Data file list
Scans = cell(1,1);
Scans(1,1) = {'3T/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
map_new = nan(96854,size(subs,2));

for s = 1:size(subs,2)

    % Find vertices for network 1, network 2, and network overlap
    subj = ft_read_cifti(sprintf('Results/Maps/Example_Subject_%d.dtseries.nii',subs(s)));
    map7 = subj.dtseries(:,7); map8 = subj.dtseries(:,8);
    other = nansum(subj.dtseries(:,setdiff(1:size(subj.dtseries,2),[7 8])),2);
    Iboth = find(map7>1 & map8>1 & other<0.1);
    I7 = find(map7>1 & map8<0.1 & other<0.1);
    I8 = find(map8>1 & map7<0.1 & other<0.1);
    Nvertices(s,:) = [length(I7) length(I8) length(Iboth)];
    clear other subj
    map_new(I7,s) = 1; map_new(I8,s) = 5; map_new(Iboth,s) = 10;

    
    for n = 1:size(Scans,1)
        fprintf('Loading scan %d for subject %d\n',n,subs(s));
        D = ft_read_cifti(sprintf('HCP_data/%d/%s',subs(s),Scans{n}));
        d = [mean(D.dtseries(I7,:))' mean(D.dtseries(I8,:))' mean(D.dtseries(Iboth,:))'];
        d = d - repmat(mean(d),size(d,1),1);
        d = d./repmat(std(d),size(d,1),1);
        data_original(n,1) = {d};
        T(n,1) = {size(d,1)};
        figure; set(gcf,'Position',[10 10 375 150],'PaperPositionMode','auto')
        plot(d(1:100,1),'k','LineWidth',2); hold on
        plot(d(1:100,2),'color',[0 0.5 0],'LineWidth',2); plot(d(1:100,3),'r','LineWidth',2); hold off
        title('Real data')
        print(gcf,'Results/HMM_simulation_explanation1','-dpng','-r300');

        % Dynamic switching hypothesis 25 TRs
        window = 25*2; N1 = d(:,1); N2 = d(:,2);
        O_switch25 = zeros(size(N1));
        O1 = nan(size(N1)); O2 = nan(size(N1));
        for tp = 1:window:size(d,1)
            O_switch25(tp:tp+window/2-1) = N1(tp:tp+window/2-1);
            O1(tp:tp+window/2-1) = N1(tp:tp+window/2-1);
            O_switch25(tp+window/2:tp+window-1) = N2(tp+window/2:tp+window-1);
            O2(tp+window/2:tp+window-1) = N2(tp+window/2:tp+window-1);
        end
        data_switch25(n) = {[N1 N2 O_switch25]};
        figure; set(gcf,'Position',[10 10 1000 600],'PaperPositionMode','auto')
        subplot(3,4,10);
        plot(O_switch25(1:100),'r'); hold on;
        plot(O1(1:100),'k','LineWidth',2); plot(O2(1:100),'color',[0 0.5 0],'LineWidth',2); hold off
        title('Switching hypothesis 25 TRs')

        % Dynamic switching hypothesis 50 TRs
        window = 50*2; N1 = d(:,1); N2 = d(:,2);
        O_switch50 = zeros(size(N1));
        O1 = nan(size(N1)); O2 = nan(size(N1));
        for tp = 1:window:size(d,1)
            O_switch50(tp:tp+window/2-1) = N1(tp:tp+window/2-1);
            O1(tp:tp+window/2-1) = N1(tp:tp+window/2-1);
            O_switch50(tp+window/2:tp+window-1) = N2(tp+window/2:tp+window-1);
            O2(tp+window/2:tp+window-1) = N2(tp+window/2:tp+window-1);
        end
        data_switch50(n) = {[N1 N2 O_switch50]};
        subplot(3,4,11);
        plot(O_switch50(1:100),'r'); hold on;
        plot(O1(1:100),'k','LineWidth',2); plot(O2(1:100),'color',[0 0.5 0],'LineWidth',2); hold off
        title('Switching hypothesis 50 TRs')

        % Dynamic switching hypothesis 10 TRs
        window = 10*2; N1 = d(:,1); N2 = d(:,2);
        O_switch10 = zeros(size(N1));
        O1 = nan(size(N1)); O2 = nan(size(N1));
        for tp = 1:window:size(d,1)
            O_switch10(tp:tp+window/2-1) = N1(tp:tp+window/2-1);
            O1(tp:tp+window/2-1) = N1(tp:tp+window/2-1);
            O_switch10(tp+window/2:tp+window-1) = N2(tp+window/2:tp+window-1);
            O2(tp+window/2:tp+window-1) = N2(tp+window/2:tp+window-1);
        end
        data_switch10(n) = {[N1 N2 O_switch10]};
        subplot(3,4,9);
        plot(O_switch10(1:100),'r'); hold on;
        plot(O1(1:100),'k','LineWidth',2); plot(O2(1:100),'color',[0 0.5 0],'LineWidth',2); hold off
        title('Switching hypothesis 10 TRs')

        % Max switching hypothesis
        [O_max,I] = max([N1 N2],[],2);
        data_max(n) = {[N1 N2 O_max]};
        O1 = nan(size(N1)); O2 = nan(size(N1));
        O1(I==1) = N1(I==1);
        O2(I==2) = N2(I==2);
        subplot(3,4,12);
        plot(O_max(1:100),'r'); hold on;
        plot(O1(1:100),'k','LineWidth',2); plot(O2(1:100),'color',[0 0.5 0],'LineWidth',2); hold off
        title('Switching hypothesis maximum')

        % Linear (additive) coupling hypothesis
        O_linear = N1+N2;
        data_linear(n) = {[N1 N2 O_linear]};
        subplot(3,4,5)
        plot(O_linear(1:100),'LineWidth',2)
        title('Coupling hypothesis: additive')

        % Non-linear (multiplicative) coupling hypothesis
        O_nonlinear = N1.*N2;
        data_nonlinear(n) = {[N1 N2 O_nonlinear]};
        subplot(3,4,6)
        plot(O_nonlinear(1:100),'LineWidth',2)
        title('Coupling hypothesis: multiplicative')

        % Spatial random mixture hypothesis
        O_mix = D.dtseries(I7(randi(length(I7),floor(length(Iboth/2)),1)),:);
        O_mix = [O_mix; D.dtseries(I8(randi(length(I8),ceil(length(Iboth/2)),1)),:)];
        O_mix = O_mix-repmat(mean(O_mix),size(O_mix,1),1);
        O_mix = O_mix./repmat(std(O_mix),size(O_mix,1),1);
        O_mix = mean(O_mix)'; O_mix = O_mix-mean(O_mix); O_mix = O_mix/std(O_mix);
        data_mix(n) = {[N1 N2 O_mix]};
        subplot(3,4,1)
        plot(O_mix(1:100),'LineWidth',2)
        title('Mixing hypothesis: random')

        % Spatial gradient hypothesis
        I7new = I7; I8new = I8;
        [~,I] = max([map7(Iboth) map8(Iboth)],[],2);
        O_gradient = zeros(size(D.dtseries,2),length(Iboth));
        for v = 1:length(Iboth)
            if I(v)==1
                [~,i] = min(I7new-I(v));
                O_gradient(:,v) = D.dtseries(I7new(i),:)';
                %I7new(i) = nan;
            elseif I(v)==2
                [~,i] = min(I8new-I(v));
                O_gradient(:,v) = D.dtseries(I8new(i),:)';
                %I8new(i) = nan;
            end
        end
        O_gradient = O_gradient-repmat(mean(O_gradient),size(O_gradient,1),1);
        O_gradient = O_gradient./repmat(std(O_gradient),size(O_gradient,1),1);
        O_gradient = mean(O_gradient,2); O_gradient = O_gradient-mean(O_gradient); O_gradient = O_gradient/std(O_gradient);
        data_gradient(n) = {[N1 N2 O_gradient]};
        subplot(3,4,2)
        plot(O_gradient(1:100),'LineWidth',2)
        title('Mixing hypothesis: gradient')

    end
end

print(gcf,'Results/HMM_simulation_explanation2','-dpng','-r300');

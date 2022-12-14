clear all; close all; clc

addpath('/Users/janinebijsterbosch/Dropbox/Matlab/FSLNets')
addpath('~/Dropbox/Matlab/hline_vline')
load('Results/HMM_output.mat','subs')

run = 'hmm'; % 'timeseries' to run part 1 or 'hmm' to run part 2

% [steve] Would be interesting to regress N1 and N2 timeseries out of O and look at residuals
% Maybe scatter plot of residual vs original (or vs N1 or N2)
% Also power spectrum plot to see if there still seems to be low frequency BOLD information
% [steve] Explore histogram/skew/spectrum of N1/ N2/ O timeseries more generally
if strcmp(run,'timeseries')
    ts_spectra_3T = nan(600,20,8); ts_spectra_7T = nan(450,20,4);
    skew_3T = nan(20,8,4); skew_7T = nan(20,4,4);
    kurt_3T = nan(20,8,4); kurt_7T = nan(20,4,4);
    for s = setdiff(1:20,[7 14])
        load(sprintf('Results/HMM_output_subject%02d_%d.mat',s,subs(s)));
        fprintf('running stats for subject %d\n',s)
        for n = 1:8
            d = data_original{n};
            Oresidual = d(:,3) - d(:,1:2)*(pinv(d(:,1:2))*d(:,3));
            spectrum = abs(fft(nets_demean(Oresidual)));
            ts_spectra_3T(:,s,n) = spectrum(1:round(size(spectrum,1)/2));
            skew_3T(s,n,:) = [skewness(d) skewness(Oresidual)];
            kurt_3T(s,n,:) = [kurtosis(d) kurtosis(Oresidual)];
        end
        for n = 9:12
            d = data_original{n};
            Oresidual = d(:,3) - d(:,1:2)*(pinv(d(:,1:2))*d(:,3));
            spectrum = abs(fft(nets_demean(Oresidual)));
            ts_spectra_7T(:,s,n-8) = spectrum(1:round(size(spectrum,1)/2));
            skew_7T(s,n-8,:) = [skewness(d) skewness(Oresidual)];
            kurt_7T(s,n-8,:) = [kurtosis(d) kurtosis(Oresidual)];
        end
    end
    skew_3T = reshape(skew_3T,20*8,4); skew_7T = reshape(skew_7T,20*4,4);
    kurt_3T = reshape(kurt_3T,20*8,4); kurt_7T = reshape(kurt_7T,20*4,4);

    % Plot spectra
    ts_spectra_3T = squeeze(nanmean(ts_spectra_3T,3));
    ts_spectra_3T = ts_spectra_3T ./ repmat(nanmax(ts_spectra_3T),size(ts_spectra_3T,1),1);
    median_spectra_3T = nanmedian(ts_spectra_3T,2);
    ts_spectra_7T = squeeze(nanmean(ts_spectra_7T,3));
    ts_spectra_7T = ts_spectra_7T ./ repmat(nanmax(ts_spectra_7T),size(ts_spectra_7T,1),1);
    median_spectra_7T = nanmedian(ts_spectra_7T,2);
    figure; set(gcf,'Position',[1000 200 750 800],'PaperPositionMode','auto')
    subplot(3,2,1); plot(ts_spectra_3T); hold on; plot(median_spectra_3T,'k','linewidth',2)
    title('Spectra of O residual timeseries in 3T data'); axis([0 600 0 1])
    subplot(3,2,2); plot(ts_spectra_7T); hold on; plot(median_spectra_7T,'k','linewidth',2)
    title('Spectra of O residual timeseries in 78T data'); axis([0 450 0 1])

    % Plot skewness and kurtosis
    subplot(3,2,3)
    swarmchart(1*ones(1,size(skew_3T,1)),skew_3T(:,1),5,[0 0 1],'filled'); hold on;
    plot([0.75 1.25],[nanmean(skew_3T(:,1)) nanmean(skew_3T(:,1))],'k'); plot([0.75 1.25],[nanmedian(skew_3T(:,1)) nanmedian(skew_3T(:,1))],'k','linewidth',2);
    swarmchart(2*ones(1,size(skew_3T,1)),skew_3T(:,2),5,[0 1 1],'filled'); hold on;
    plot([1.75 2.25],[nanmean(skew_3T(:,2)) nanmean(skew_3T(:,2))],'k'); plot([1.75 2.25],[nanmedian(skew_3T(:,2)) nanmedian(skew_3T(:,2))],'k','linewidth',2);
    swarmchart(3*ones(1,size(skew_3T,1)),skew_3T(:,3),5,[0 0.4470 0.7410],'filled'); hold on;
    plot([2.75 3.25],[nanmean(skew_3T(:,3)) nanmean(skew_3T(:,3))],'k'); plot([2.75 3.25],[nanmedian(skew_3T(:,3)) nanmedian(skew_3T(:,3))],'k','linewidth',2);
    swarmchart(4*ones(1,size(skew_3T,1)),skew_3T(:,4),5,[0.4940 0.1840 0.5560],'filled'); hold on;
    plot([3.75 4.25],[nanmean(skew_3T(:,4)) nanmean(skew_3T(:,4))],'k'); plot([3.75 4.25],[nanmedian(skew_3T(:,4)) nanmedian(skew_3T(:,4))],'k','linewidth',2);
    axis([0.5 4.5 -0.8 0.8]); hline(0,'-k')
    set(gca,'xtick',1:1:4,'xticklabel',{'N1','N2','O','O residuals'});
    ylabel('Skewness'); title('Timeseries skewness in 3T data')
    subplot(3,2,4)
    swarmchart(1*ones(1,size(skew_7T,1)),skew_7T(:,1),5,[0 0 1],'filled'); hold on;
    plot([0.75 1.25],[nanmean(skew_7T(:,1)) nanmean(skew_7T(:,1))],'k'); plot([0.75 1.25],[nanmedian(skew_7T(:,1)) nanmedian(skew_7T(:,1))],'k','linewidth',2);
    swarmchart(2*ones(1,size(skew_7T,1)),skew_7T(:,2),5,[0 1 1],'filled'); hold on;
    plot([1.75 2.25],[nanmean(skew_7T(:,2)) nanmean(skew_7T(:,2))],'k'); plot([1.75 2.25],[nanmedian(skew_7T(:,2)) nanmedian(skew_7T(:,2))],'k','linewidth',2);
    swarmchart(3*ones(1,size(skew_7T,1)),skew_7T(:,3),5,[0 0.4470 0.7410],'filled'); hold on;
    plot([2.75 3.25],[nanmean(skew_7T(:,3)) nanmean(skew_7T(:,3))],'k'); plot([2.75 3.25],[nanmedian(skew_7T(:,3)) nanmedian(skew_7T(:,3))],'k','linewidth',2);
    swarmchart(4*ones(1,size(skew_7T,1)),skew_7T(:,4),5,[0.4940 0.1840 0.5560],'filled'); hold on;
    plot([3.75 4.25],[nanmean(skew_7T(:,4)) nanmean(skew_7T(:,4))],'k'); plot([3.75 4.25],[nanmedian(skew_7T(:,4)) nanmedian(skew_7T(:,4))],'k','linewidth',2);
    axis([0.5 4.5 -0.8 0.8]); hline(0,'-k')
    set(gca,'xtick',1:1:4,'xticklabel',{'N1','N2','O','O residuals'});
    ylabel('Skewness'); title('Timeseries skewness in 7T data')
    subplot(3,2,5)
    swarmchart(1*ones(1,size(kurt_3T,1)),kurt_3T(:,1),5,[0 0 1],'filled'); hold on;
    plot([0.75 1.25],[nanmean(kurt_3T(:,1)) nanmean(kurt_3T(:,1))],'k'); plot([0.75 1.25],[nanmedian(kurt_3T(:,1)) nanmedian(kurt_3T(:,1))],'k','linewidth',2);
    swarmchart(2*ones(1,size(kurt_3T,1)),kurt_3T(:,2),5,[0 1 1],'filled'); hold on;
    plot([1.75 2.25],[nanmean(kurt_3T(:,2)) nanmean(kurt_3T(:,2))],'k'); plot([1.75 2.25],[nanmedian(kurt_3T(:,2)) nanmedian(kurt_3T(:,2))],'k','linewidth',2);
    swarmchart(3*ones(1,size(kurt_3T,1)),kurt_3T(:,3),5,[0 0.4470 0.7410],'filled'); hold on;
    plot([2.75 3.25],[nanmean(kurt_3T(:,3)) nanmean(kurt_3T(:,3))],'k'); plot([2.75 3.25],[nanmedian(kurt_3T(:,3)) nanmedian(kurt_3T(:,3))],'k','linewidth',2);
    swarmchart(4*ones(1,size(kurt_3T,1)),kurt_3T(:,4),5,[0.4940 0.1840 0.5560],'filled'); hold on;
    plot([3.75 4.25],[nanmean(kurt_3T(:,4)) nanmean(kurt_3T(:,4))],'k'); plot([3.75 4.25],[nanmedian(kurt_3T(:,4)) nanmedian(kurt_3T(:,4))],'k','linewidth',2);
    axis([0.5 4.5 2 5])
    set(gca,'xtick',1:1:4,'xticklabel',{'N1','N2','O','O residuals'});
    ylabel('Kurtosis'); title('Timeseries kurtosis in 3T data')
    subplot(3,2,6)
    swarmchart(1*ones(1,size(kurt_7T,1)),kurt_7T(:,1),5,[0 0 1],'filled'); hold on;
    plot([0.75 1.25],[nanmean(kurt_7T(:,1)) nanmean(kurt_7T(:,1))],'k'); plot([0.75 1.25],[nanmedian(kurt_7T(:,1)) nanmedian(kurt_7T(:,1))],'k','linewidth',2);
    swarmchart(2*ones(1,size(kurt_7T,1)),kurt_7T(:,2),5,[0 1 1],'filled'); hold on;
    plot([1.75 2.25],[nanmean(kurt_7T(:,2)) nanmean(kurt_7T(:,2))],'k'); plot([1.75 2.25],[nanmedian(kurt_7T(:,2)) nanmedian(kurt_7T(:,2))],'k','linewidth',2);
    swarmchart(3*ones(1,size(kurt_7T,1)),kurt_7T(:,3),5,[0 0.4470 0.7410],'filled'); hold on;
    plot([2.75 3.25],[nanmean(kurt_7T(:,3)) nanmean(kurt_7T(:,3))],'k'); plot([2.75 3.25],[nanmedian(kurt_7T(:,3)) nanmedian(kurt_7T(:,3))],'k','linewidth',2);
    swarmchart(4*ones(1,size(kurt_7T,1)),kurt_7T(:,4),5,[0.4940 0.1840 0.5560],'filled'); hold on;
    plot([3.75 4.25],[nanmean(kurt_7T(:,4)) nanmean(kurt_7T(:,4))],'k'); plot([3.75 4.25],[nanmedian(kurt_7T(:,4)) nanmedian(kurt_7T(:,4))],'k','linewidth',2);
    axis([0.5 4.5 2 5])
    set(gca,'xtick',1:1:4,'xticklabel',{'N1','N2','O','O residuals'});
    ylabel('Kurtosis'); title('Timeseries kurtosis in 7T data')
    print(gcf,'Results/October_feedback_skewness','-dpng','-r300');

    %% HMM summary measures
elseif strcmp(run,'hmm')
    load('Results/HMM_output.mat','StateTransitions_original','LifeTimes_original','LifeTimes_linear','StateTransitions_linear');
    LifeTimes_state1 = [];
    LifeTimes_state2 = [];
    LifeTimes_state3 = [];
    StateTransitions = nan(3,3,20);
    LifeTimes_state1_linear = [];
    LifeTimes_state2_linear = [];
    LifeTimes_state3_linear = [];
    StateTransitions2 = nan(3,3,20);
    for s = setdiff(1:length(subs),[7 14])
        LT = LifeTimes_original(s,:);
        LifeTimes_state1 = [LifeTimes_state1;LT{1}'];
        LifeTimes_state2 = [LifeTimes_state2;LT{2}'];
        LifeTimes_state3 = [LifeTimes_state3;LT{3}'];
        StateTransitions(:,:,s) = StateTransitions_original{s};

        LT = LifeTimes_linear(s,:);
        LifeTimes_state1_linear = [LifeTimes_state1_linear;LT{1}'];
        LifeTimes_state2_linear = [LifeTimes_state2_linear;LT{2}'];
        LifeTimes_state3_linear = [LifeTimes_state3_linear;LT{3}'];
        StateTransitions2(:,:,s) = StateTransitions_linear{s};
    end

    figure; set(gcf,'Position',[600 600 1000 400],'PaperPositionMode','auto')
    subplot(2,2,1); histogram(LifeTimes_state1); hold on; histogram(LifeTimes_state2); histogram(LifeTimes_state3);
    legend({'State1','State2','State3'}); title('HMM state lifetimes original'); axis([0 100 0 1000])
    subplot(2,2,2); histogram(LifeTimes_state1_linear); hold on; histogram(LifeTimes_state2_linear); histogram(LifeTimes_state3_linear);
    legend({'State1','State2','State3'}); title('HMM state lifetimes linear simulated data'); axis([0 100 0 1000])
    subplot(2,2,3); imagesc(nanmean(StateTransitions,3),[0 1]); colorbar; title('HMM state transitions original')
    subplot(2,2,4); imagesc(nanmean(StateTransitions2,3),[0 1]); colorbar; title('HMM state transitions linear simulated data')
    
    print(gcf,'Results/October_feedback_hmm','-dpng','-r300');

end

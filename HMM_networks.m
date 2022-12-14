clear all; close all; clc

% Set paths
addpath(genpath('~/Dropbox/Matlab/HMM-MAR-master/'))
addpath(genpath('~/Dropbox/Matlab/fieldtrip-20190819/'),'-END')
addpath('~/Dropbox/Matlab/hline_vline')

% Set options for HMM
options.Fs = 1/0.72;
options.order = 0;
options.covtype = 'uniquefull';
repetitions = 5;
maxk = 4;
RunHMM = 1;

% Subject list
subs = [105923 114823 125525 130518 137128 144226 146129 158035 169343 177746 185442 192439 195041 200614 204521 562345 601127 783462 859671 861456];

% Data file list
Scans = cell(12,1);
Scans(1,1) = {'3T/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
Scans(2,1) = {'3T/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
Scans(3,1) = {'3T/rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
Scans(4,1) = {'3T/rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
Scans(5,1) = {'retest/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
Scans(6,1) = {'retest/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
Scans(7,1) = {'retest/rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
Scans(8,1) = {'retest/rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
Scans(9,1) = {'7T/rfMRI_REST1_7T_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
Scans(10,1) = {'7T/rfMRI_REST2_7T_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
Scans(11,1) = {'7T/rfMRI_REST3_7T_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
Scans(12,1) = {'7T/rfMRI_REST4_7T_Atlas_MSMAll_hp2000_clean.dtseries.nii'};

Nvertices = zeros(length(subs),3);
TScorr = nan(size(subs,2)*12,8);
Means_original_all = cell(length(subs),maxk-1);
Means_switch50_all = cell(length(subs),maxk-1);
Means_switch25_all = cell(length(subs),maxk-1);
Means_switch10_all = cell(length(subs),maxk-1);
Means_linear_all = cell(length(subs),maxk-1);
Means_max_all = cell(length(subs),maxk-1);
Means_nonlinear_all = cell(length(subs),maxk-1);
Means_mix_all = cell(length(subs),maxk-1);
Means_gradient_all = cell(length(subs),maxk-1);
GammaSim_original = nan(size(subs,2)*(repetitions*(repetitions-1)/2),maxk-1);
GammaSim_switch25 = nan(size(subs,2)*(repetitions*(repetitions-1)/2),maxk-1); MeansSim_switch25 = nan(size(subs,2)*repetitions*repetitions,maxk-1);
GammaSim_switch50 = nan(size(subs,2)*(repetitions*(repetitions-1)/2),maxk-1); MeansSim_switch50 = nan(size(subs,2)*repetitions*repetitions,maxk-1);
GammaSim_switch10 = nan(size(subs,2)*(repetitions*(repetitions-1)/2),maxk-1); MeansSim_switch10 = nan(size(subs,2)*repetitions*repetitions,maxk-1);
GammaSim_linear = nan(size(subs,2)*(repetitions*(repetitions-1)/2),maxk-1); MeansSim_linear = nan(size(subs,2)*repetitions*repetitions,maxk-1);
GammaSim_nonlinear = nan(size(subs,2)*(repetitions*(repetitions-1)/2),maxk-1); MeansSim_nonlinear = nan(size(subs,2)*repetitions*repetitions,maxk-1);
GammaSim_max = nan(size(subs,2)*(repetitions*(repetitions-1)/2),maxk-1); MeansSim_max = nan(size(subs,2)*repetitions*repetitions,maxk-1);
GammaSim_gradient = nan(size(subs,2)*(repetitions*(repetitions-1)/2),maxk-1); MeansSim_gradient = nan(size(subs,2)*repetitions*repetitions,maxk-1);
GammaSim_mix = nan(size(subs,2)*(repetitions*(repetitions-1)/2),maxk-1); MeansSim_mix = nan(size(subs,2)*repetitions*repetitions,maxk-1);
map_new = nan(96854,size(subs,2));
LifeTimes_original = cell(length(subs),3);
StateTransitions_original = cell(length(subs),1);
LifeTimes_linear = cell(length(subs),3);
StateTransitions_linear = cell(length(subs),1);

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

    if RunHMM
        if Iboth

            % Prep data for HMM
            data_original = cell(size(Scans,1),1);
            data_switch25 = cell(size(Scans,1),1);
            data_switch50 = cell(size(Scans,1),1);
            data_switch10 = cell(size(Scans,1),1);
            data_linear = cell(size(Scans,1),1);
            data_nonlinear = cell(size(Scans,1),1);
            data_max = cell(size(Scans,1),1);
            data_gradient = cell(size(Scans,1),1);
            data_mix = cell(size(Scans,1),1);

            T = cell(size(Scans,1),1);
            Overlap_correlations = [];
            for n = 1:size(Scans,1)
                fprintf('Loading scan %d for subject %d\n',n,subs(s));
                D = ft_read_cifti(sprintf('HCP_data/%d/%s',subs(s),Scans{n}));
                d = [mean(D.dtseries(I7,:))' mean(D.dtseries(I8,:))' mean(D.dtseries(Iboth,:))'];
                d = d - repmat(mean(d),size(d,1),1);
                d = d./repmat(std(d),size(d,1),1);
                data_original(n,1) = {d};
                T(n,1) = {size(d,1)};

                % Dynamic switching hypothesis 25 TRs
                window = 25*2; N1 = d(:,1); N2 = d(:,2);
                O_switch25 = zeros(size(N1));
                for tp = 1:window:size(d,1)
                    O_switch25(tp:tp+window/2-1) = N1(tp:tp+window/2-1);
                    O_switch25(tp+window/2:tp+window-1) = N2(tp+window/2:tp+window-1);
                end
                data_switch25(n) = {[N1 N2 O_switch25]};

                % Dynamic switching hypothesis 50 TRs
                window = 50*2; N1 = d(:,1); N2 = d(:,2);
                O_switch50 = zeros(size(N1));
                for tp = 1:window:size(d,1)
                    O_switch50(tp:tp+window/2-1) = N1(tp:tp+window/2-1);
                    O_switch50(tp+window/2:tp+window-1) = N2(tp+window/2:tp+window-1);
                end
                data_switch50(n) = {[N1 N2 O_switch50]};

                % Dynamic switching hypothesis 10 TRs
                window = 10*2; N1 = d(:,1); N2 = d(:,2);
                O_switch10 = zeros(size(N1));
                for tp = 1:window:size(d,1)
                    O_switch10(tp:tp+window/2-1) = N1(tp:tp+window/2-1);
                    O_switch10(tp+window/2:tp+window-1) = N2(tp+window/2:tp+window-1);
                end
                data_switch10(n) = {[N1 N2 O_switch10]};

                % Max switching hypothesis
                O_max = max([N1 N2],[],2);
                data_max(n) = {[N1 N2 O_max]};

                % Linear (additive) coupling hypothesis
                O_linear = N1+N2;
                data_linear(n) = {[N1 N2 O_linear]};

                % Non-linear (multiplicative) coupling hypothesis
                O_nonlinear = N1.*N2;
                data_nonlinear(n) = {[N1 N2 O_nonlinear]};

                % Spatial random mixture hypothesis
                O_mix = D.dtseries(I7(randi(length(I7),floor(length(Iboth/2)),1)),:);
                O_mix = [O_mix; D.dtseries(I8(randi(length(I8),ceil(length(Iboth/2)),1)),:)];
                O_mix = O_mix-repmat(mean(O_mix),size(O_mix,1),1);
                O_mix = O_mix./repmat(std(O_mix),size(O_mix,1),1);
                O_mix = mean(O_mix)'; O_mix = O_mix-mean(O_mix); O_mix = O_mix/std(O_mix);
                data_mix(n) = {[N1 N2 O_mix]};

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

                r = corr([d(:,3) O_switch50 O_switch25 O_switch10 O_max O_linear O_nonlinear O_mix O_gradient]);
                Overlap_correlations = [Overlap_correlations; r(1,2:end)];
            end
            clear n I1 I2 Iboth
            TScorr((s-1)*12+1:s*12,:) = Overlap_correlations;

            % Run HMM
            for k = 2:maxk
                options.K = k;
                Gamma_originals = cell(repetitions,1); Means_original = cell(repetitions,1);
                Gamma_switch25s = cell(repetitions,1); Means_switch25 = cell(repetitions,1);
                Gamma_switch50s = cell(repetitions,1); Means_switch50 = cell(repetitions,1);
                Gamma_switch10s = cell(repetitions,1); Means_switch10 = cell(repetitions,1);
                Gamma_linears = cell(repetitions,1); Means_linear = cell(repetitions,1);
                Gamma_nonlinears = cell(repetitions,1); Means_nonlinear = cell(repetitions,1);
                Gamma_maxs = cell(repetitions,1); Means_max = cell(repetitions,1);
                Gamma_gradients = cell(repetitions,1); Means_gradient = cell(repetitions,1);
                Gamma_mixs = cell(repetitions,1); Means_mix = cell(repetitions,1);
                for reps = 1:repetitions
                    [hmm_original, Gamma_original,~,~,~,~,~] = hmmmar(data_original,T,options);
                    Gamma_originals(reps) = {Gamma_original}; Means_original(reps) = {getMean(hmm_original)};
                    if k==3 && reps==1
                        LifeTimes_original(s,:) = getStateLifeTimes(Gamma_original,T,options);
                        StateTransitions_original(s) = {hmm_original.P};
                    end

                    [hmm_switch25, Gamma_switch25,~,~,~,~,~] = hmmmar(data_switch25,T,options);
                    Gamma_switch25s(reps) = {Gamma_switch25}; Means_switch25(reps) = {getMean(hmm_switch25)};

                    [hmm_switch50, Gamma_switch50,~,~,~,~,~] = hmmmar(data_switch50,T,options);
                    Gamma_switch50s(reps) = {Gamma_switch50}; Means_switch50(reps) = {getMean(hmm_switch50)};

                    [hmm_switch10, Gamma_switch10,~,~,~,~,~] = hmmmar(data_switch10,T,options);
                    Gamma_switch10s(reps) = {Gamma_switch10}; Means_switch10(reps) = {getMean(hmm_switch10)};

                    [hmm_linear, Gamma_linear,~,~,~,~,~] = hmmmar(data_linear,T,options);
                    Gamma_linears(reps) = {Gamma_linear}; Means_linear(reps) = {getMean(hmm_linear)};
                    if k==3 && reps==1
                        LifeTimes_linear(s,:) = getStateLifeTimes(Gamma_linear,T,options);
                        StateTransitions_linear(s) = {hmm_linear.P};
                    end

                    [hmm_nonlinear, Gamma_nonlinear,~,~,~,~,~] = hmmmar(data_nonlinear,T,options);
                    Gamma_nonlinears(reps) = {Gamma_nonlinear}; Means_nonlinear(reps) = {getMean(hmm_nonlinear)};

                    [hmm_max, Gamma_max,~,~,~,~,~] = hmmmar(data_max,T,options);
                    Gamma_maxs(reps) = {Gamma_max}; Means_max(reps) = {getMean(hmm_max)};

                    [hmm_gradient, Gamma_gradient,~,~,~,~,~] = hmmmar(data_gradient,T,options);
                    Gamma_gradients(reps) = {Gamma_gradient}; Means_gradient(reps) = {getMean(hmm_gradient)};

                    [hmm_mix, Gamma_mix,~,~,~,~,~] = hmmmar(data_mix,T,options);
                    Gamma_mixs(reps) = {Gamma_mix}; Means_mix(reps) = {getMean(hmm_mix)};
                end
                i=1;
                for rep1 = 1:repetitions
                    for rep2 = 1:repetitions
                        if rep2>rep1
                            GammaSim_original((s-1)*(repetitions*(repetitions-1)/2)+i,k-1) = getGammaSimilarity(Gamma_originals{rep1},Gamma_originals{rep2});
                            GammaSim_switch25((s-1)*(repetitions*(repetitions-1)/2)+i,k-1) = getGammaSimilarity(Gamma_switch25s{rep1},Gamma_switch25s{rep2});
                            GammaSim_switch50((s-1)*(repetitions*(repetitions-1)/2)+i,k-1) = getGammaSimilarity(Gamma_switch50s{rep1},Gamma_switch50s{rep2});
                            GammaSim_switch10((s-1)*(repetitions*(repetitions-1)/2)+i,k-1) = getGammaSimilarity(Gamma_switch10s{rep1},Gamma_switch10s{rep2});
                            GammaSim_linear((s-1)*(repetitions*(repetitions-1)/2)+i,k-1) = getGammaSimilarity(Gamma_linears{rep1},Gamma_linears{rep2});
                            GammaSim_nonlinear((s-1)*(repetitions*(repetitions-1)/2)+i,k-1) = getGammaSimilarity(Gamma_nonlinears{rep1},Gamma_nonlinears{rep2});
                            GammaSim_max((s-1)*(repetitions*(repetitions-1)/2)+i,k-1) = getGammaSimilarity(Gamma_maxs{rep1},Gamma_maxs{rep2});
                            GammaSim_gradient((s-1)*(repetitions*(repetitions-1)/2)+i,k-1) = getGammaSimilarity(Gamma_gradients{rep1},Gamma_gradients{rep2});
                            GammaSim_mix((s-1)*(repetitions*(repetitions-1)/2)+i,k-1) = getGammaSimilarity(Gamma_mixs{rep1},Gamma_mixs{rep2});
                            i = i+1;
                        end
                    end
                end

                % Align means for all repetitions
                for rep2 = 2:repetitions
                    [~,assign] = matrixICC(Means_original{1},Means_original{rep2},1); Means_original(rep2) = {Means_original{rep2}(:,assign)};
                    [~,assign] = matrixICC(Means_switch50{1},Means_switch50{rep2},1); Means_switch50(rep2) = {Means_switch50{rep2}(:,assign)};
                    [~,assign] = matrixICC(Means_switch25{1},Means_switch25{rep2},1); Means_switch25(rep2) = {Means_switch25{rep2}(:,assign)};
                    [~,assign] = matrixICC(Means_switch10{1},Means_switch10{rep2},1); Means_switch10(rep2) = {Means_switch10{rep2}(:,assign)};
                    [~,assign] = matrixICC(Means_max{1},Means_max{rep2},1); Means_max(rep2) = {Means_max{rep2}(:,assign)};
                    [~,assign] = matrixICC(Means_linear{1},Means_linear{rep2},1); Means_linear(rep2) = {Means_linear{rep2}(:,assign)};
                    [~,assign] = matrixICC(Means_nonlinear{1},Means_nonlinear{rep2},1); Means_nonlinear(rep2) = {Means_nonlinear{rep2}(:,assign)};
                    [~,assign] = matrixICC(Means_mix{1},Means_mix{rep2},1); Means_mix(rep2) = {Means_mix{rep2}(:,assign)};
                    [~,assign] = matrixICC(Means_gradient{1},Means_gradient{rep2},1); Means_gradient(rep2) = {Means_gradient{rep2}(:,assign)};
                end

                % Calculate mean of state means across repetitions
                M_original = []; M_switch50 = []; M_switch25 = []; M_switch10 = []; M_max = []; M_linear = []; M_nonlinear = []; M_mix = []; M_gradient = [];
                for rep1 = 1:repetitions
                    M_original = cat(3,M_original,Means_original{rep1});
                    M_switch50 = cat(3,M_switch50,Means_switch50{rep1});
                    M_switch25 = cat(3,M_switch25,Means_switch25{rep1});
                    M_switch10 = cat(3,M_switch10,Means_switch10{rep1});
                    M_max = cat(3,M_max,Means_max{rep1});
                    M_linear = cat(3,M_linear,Means_linear{rep1});
                    M_nonlinear = cat(3,M_nonlinear,Means_nonlinear{rep1});
                    M_mix = cat(3,M_mix,Means_mix{rep1});
                    M_gradient = cat(3,M_gradient,Means_gradient{rep1});
                end
                Means_original_all{s,k-1} = mean(M_original,3);
                Means_switch50_all{s,k-1} = mean(M_switch50,3);
                Means_switch25_all{s,k-1} = mean(M_switch25,3);
                Means_switch10_all{s,k-1} = mean(M_switch10,3);
                Means_max_all{s,k-1} = mean(M_max,3);
                Means_linear_all{s,k-1} = mean(M_linear,3);
                Means_nonlinear_all{s,k-1} = mean(M_nonlinear,3);
                Means_mix_all{s,k-1} = mean(M_mix,3);
                Means_gradient_all{s,k-1} = mean(M_gradient,3);

                % Align means with original
                [~,assign] = matrixICC(Means_original_all{s,k-1},Means_switch25_all{s,k-1},1);
                Means_switch25_all{s,k-1} = Means_switch25_all{s,k-1}(:,assign);
                for rep1 = 1:repetitions;  Means_switch25{rep1} = Means_switch25{rep1}(:,assign); end

                [~,assign] = matrixICC(Means_original_all{s,k-1},Means_switch10_all{s,k-1},1);
                Means_switch10_all{s,k-1} = Means_switch10_all{s,k-1}(:,assign);
                for rep1 = 1:repetitions;  Means_switch10{rep1} = Means_switch10{rep1}(:,assign); end

                [~,assign] = matrixICC(Means_original_all{s,k-1},Means_max_all{s,k-1},1);
                Means_max_all{s,k-1} = Means_max_all{s,k-1}(:,assign);
                for rep1 = 1:repetitions;  Means_max{rep1} = Means_max{rep1}(:,assign); end

                [~,assign] = matrixICC(Means_original_all{s,k-1},Means_switch50_all{s,k-1},1);
                Means_switch50_all{s,k-1} = Means_switch50_all{s,k-1}(:,assign);
                for rep1 = 1:repetitions;  Means_switch50{rep1} = Means_switch50{rep1}(:,assign); end

                [~,assign] = matrixICC(Means_original_all{s,k-1},Means_linear_all{s,k-1},1);
                Means_linear_all{s,k-1} = Means_linear_all{s,k-1}(:,assign);
                for rep1 = 1:repetitions;  Means_linear{rep1} = Means_linear{rep1}(:,assign); end

                [~,assign] = matrixICC(Means_original_all{s,k-1},Means_nonlinear_all{s,k-1},1);
                Means_nonlinear_all{s,k-1} = Means_nonlinear_all{s,k-1}(:,assign);
                for rep1 = 1:repetitions;  Means_nonlinear{rep1} = Means_nonlinear{rep1}(:,assign); end

                [~,assign] = matrixICC(Means_original_all{s,k-1},Means_mix_all{s,k-1},1);
                Means_mix_all{s,k-1} = Means_mix_all{s,k-1}(:,assign);
                for rep1 = 1:repetitions;  Means_mix{rep1} = Means_mix{rep1}(:,assign); end

                [~,assign] = matrixICC(Means_original_all{s,k-1},Means_gradient_all{s,k-1},1);
                Means_gradient_all{s,k-1} = Means_gradient_all{s,k-1}(:,assign);
                for rep1 = 1:repetitions;  Means_gradient{rep1} = Means_gradient{rep1}(:,assign); end

                % Correlate means of all hypotheses with original
                i=1;
                for rep1 = 1:repetitions
                    for rep2 = 1:repetitions
                        [r,~] = matrixICC(Means_original{rep1},Means_switch50{rep2},1); r = r(eye(options.K)==1);
                        MeansSim_switch50((s-1)*(repetitions*repetitions)+i,k-1) = mean(r);

                        [r,~] = matrixICC(Means_original{rep1},Means_switch25{rep2},1); r = r(eye(options.K)==1);
                        MeansSim_switch25((s-1)*(repetitions*repetitions)+i,k-1) = mean(r);

                        [r,~] = matrixICC(Means_original{rep1},Means_switch10{rep2},1); r = r(eye(options.K)==1);
                        MeansSim_switch10((s-1)*(repetitions*repetitions)+i,k-1) = mean(r);

                        [r,~] = matrixICC(Means_original{rep1},Means_max{rep2},1); r = r(eye(options.K)==1);
                        MeansSim_max((s-1)*(repetitions*repetitions)+i,k-1) = mean(r);

                        [r,~] = matrixICC(Means_original{rep1},Means_linear{rep2},1); r = r(eye(options.K)==1);
                        MeansSim_linear((s-1)*(repetitions*repetitions)+i,k-1) = mean(r);

                        [r,~] = matrixICC(Means_original{rep1},Means_nonlinear{rep2},1); r = r(eye(options.K)==1);
                        MeansSim_nonlinear((s-1)*(repetitions*repetitions)+i,k-1) = mean(r);

                        [r,~] = matrixICC(Means_original{rep1},Means_mix{rep2},1); r = r(eye(options.K)==1);
                        MeansSim_mix((s-1)*(repetitions*repetitions)+i,k-1) = mean(r);

                        [r,~] = matrixICC(Means_original{rep1},Means_gradient{rep2},1); r = r(eye(options.K)==1);
                        MeansSim_gradient((s-1)*(repetitions*repetitions)+i,k-1) = mean(r);

                        i = i+1;
                    end
                end
            end
            save(sprintf('Results/HMM_output_subject%02d_%d.mat',s,subs(s)));
        end
    end   
end

cifti = ft_read_cifti(sprintf('Results/Maps/Example_Subject_%d.dtseries.nii',subs(s)));
cifti.dtseries = map_new;
cifti.hdr.dim(7) = size(subs,2);
cifti.time = 1:size(subs,2);
ft_write_cifti('Results/Maps/Two_network_overlap',cifti,'parameter','dtseries')
clear cifti

figure; set(gcf,'Position',[1 60 1300 450],'PaperPositionMode','auto')
subplot(2,5,1);
for k = 2:maxk
    swarmchart(k*ones(1,size(GammaSim_original,1)),GammaSim_original(:,k-1),5,'k','filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(GammaSim_original(:,k-1)) nanmean(GammaSim_original(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(GammaSim_original(:,k-1)) nanmedian(GammaSim_original(:,k-1))],'k','linewidth',2)
end
title('A. Original data'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Gamma similarity'); xlabel('Number of states')
subplot(2,5,3);
for k = 2:maxk
    swarmchart(k*ones(1,size(GammaSim_original,1)),GammaSim_switch25(:,k-1),5,[0 1 1],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(GammaSim_switch25(:,k-1)) nanmean(GammaSim_switch25(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(GammaSim_switch25(:,k-1)) nanmedian(GammaSim_switch25(:,k-1))],'k','linewidth',2)
end
title('C. Switch 25 TRs'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Gamma similarity'); xlabel('Number of states')
subplot(2,5,2);
for k = 2:maxk
    swarmchart(k*ones(1,size(GammaSim_original,1)),GammaSim_switch50(:,k-1),5,[0 0 1],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(GammaSim_switch50(:,k-1)) nanmean(GammaSim_switch50(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(GammaSim_switch50(:,k-1)) nanmedian(GammaSim_switch50(:,k-1))],'k','linewidth',2)
end
title('B. Switch 50 TRs'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Gamma similarity'); xlabel('Number of states')
subplot(2,5,4);
for k = 2:maxk
    swarmchart(k*ones(1,size(GammaSim_original,1)),GammaSim_switch10(:,k-1),5,[0 0.4470 0.7410],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(GammaSim_switch10(:,k-1)) nanmean(GammaSim_switch10(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(GammaSim_switch10(:,k-1)) nanmedian(GammaSim_switch10(:,k-1))],'k','linewidth',2)
end
title('D. Switch 10 TRs'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Gamma similarity'); xlabel('Number of states')
subplot(2,5,7);
for k = 2:maxk
    swarmchart(k*ones(1,size(GammaSim_original,1)),GammaSim_linear(:,k-1),5,[0 1 0],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(GammaSim_linear(:,k-1)) nanmean(GammaSim_linear(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(GammaSim_linear(:,k-1)) nanmedian(GammaSim_linear(:,k-1))],'k','linewidth',2)
end
title('F. Linear additive coupling'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Gamma similarity'); xlabel('Number of states')
subplot(2,5,8);
for k = 2:maxk
    swarmchart(k*ones(1,size(GammaSim_original,1)),GammaSim_nonlinear(:,k-1),5,[0.4660 0.6740 0.1880],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(GammaSim_nonlinear(:,k-1)) nanmean(GammaSim_nonlinear(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(GammaSim_nonlinear(:,k-1)) nanmedian(GammaSim_nonlinear(:,k-1))],'k','linewidth',2)
end
title('G. Nonlinear multiplicative coupling'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Gamma similarity'); xlabel('Number of states')
subplot(2,5,5);
for k = 2:maxk
    swarmchart(k*ones(1,size(GammaSim_original,1)),GammaSim_max(:,k-1),5,[0.4940 0.1840 0.5560],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(GammaSim_max(:,k-1)) nanmean(GammaSim_max(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(GammaSim_max(:,k-1)) nanmedian(GammaSim_max(:,k-1))],'k','linewidth',2)
end
title('E. Max switching'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Gamma similarity'); xlabel('Number of states')
subplot(2,5,9);
for k = 2:maxk
    swarmchart(k*ones(1,size(GammaSim_original,1)),GammaSim_mix(:,k-1),5,[0.8500 0.3250 0.0980],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(GammaSim_mix(:,k-1)) nanmean(GammaSim_mix(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(GammaSim_mix(:,k-1)) nanmedian(GammaSim_mix(:,k-1))],'k','linewidth',2)
end
title('H. Spatial random mixture'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Gamma similarity'); xlabel('Number of states')
subplot(2,5,10);
for k = 2:maxk
    swarmchart(k*ones(1,size(GammaSim_original,1)),GammaSim_gradient(:,k-1),5,[0.6350 0.0780 0.1840],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(GammaSim_gradient(:,k-1)) nanmean(GammaSim_gradient(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(GammaSim_gradient(:,k-1)) nanmedian(GammaSim_gradient(:,k-1))],'k','linewidth',2)
end
title('I. Spatial gradient'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Gamma similarity'); xlabel('Number of states')
print(gcf,'Results/HMM_Gamma_similarity','-dpng','-r300');

figure; set(gcf,'Position',[620 60 1000 450],'PaperPositionMode','auto')
subplot(2,4,2);
for k = 2:maxk
    swarmchart(k*ones(1,size(MeansSim_switch25,1)),MeansSim_switch25(:,k-1),5,[0 1 1],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(MeansSim_switch25(:,k-1)) nanmean(MeansSim_switch25(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(MeansSim_switch25(:,k-1)) nanmedian(MeansSim_switch25(:,k-1))],'k','linewidth',2)
end
title('B. Switch 25 TRs'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Absolute agreement with original'); xlabel('Number of states')
subplot(2,4,1);
for k = 2:maxk
    swarmchart(k*ones(1,size(MeansSim_switch25,1)),MeansSim_switch50(:,k-1),5,[0 0.4470 0.7410],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(MeansSim_switch50(:,k-1)) nanmean(MeansSim_switch50(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(MeansSim_switch50(:,k-1)) nanmedian(MeansSim_switch50(:,k-1))],'k','linewidth',2)
end
title('A. Switch 50 TRs'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Absolute agreement with original'); xlabel('Number of states')
subplot(2,4,3);
for k = 2:maxk
    swarmchart(k*ones(1,size(MeansSim_switch25,1)),MeansSim_switch10(:,k-1),5,[0 0 1],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(MeansSim_switch10(:,k-1)) nanmean(MeansSim_switch10(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(MeansSim_switch10(:,k-1)) nanmedian(MeansSim_switch10(:,k-1))],'k','linewidth',2)
end
title('C. Switch 10 TRs'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Absolute agreement with original'); xlabel('Number of states')
subplot(2,4,5);
for k = 2:maxk
    swarmchart(k*ones(1,size(MeansSim_switch25,1)),MeansSim_linear(:,k-1),5,[0 1 0],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(MeansSim_linear(:,k-1)) nanmean(MeansSim_linear(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(MeansSim_linear(:,k-1)) nanmedian(MeansSim_linear(:,k-1))],'k','linewidth',2)
end
title('E. Linear additive coupling'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Absolute agreement with original'); xlabel('Number of states')
subplot(2,4,6);
for k = 2:maxk
    swarmchart(k*ones(1,size(MeansSim_switch25,1)),MeansSim_nonlinear(:,k-1),5,[0.4660 0.6740 0.1880],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(MeansSim_nonlinear(:,k-1)) nanmean(MeansSim_nonlinear(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(MeansSim_nonlinear(:,k-1)) nanmedian(MeansSim_nonlinear(:,k-1))],'k','linewidth',2)
end
title('F. Nonlinear multiplicative coupling'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Absolute agreement with original'); xlabel('Number of states')
subplot(2,4,4);
for k = 2:maxk
    swarmchart(k*ones(1,size(MeansSim_switch25,1)),MeansSim_max(:,k-1),5,[0.4940 0.1840 0.5560],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(MeansSim_max(:,k-1)) nanmean(MeansSim_max(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(MeansSim_max(:,k-1)) nanmedian(MeansSim_max(:,k-1))],'k','linewidth',2)
end
title('D. Max switching'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Absolute agreement with original'); xlabel('Number of states')
subplot(2,4,7);
for k = 2:maxk
    swarmchart(k*ones(1,size(MeansSim_switch25,1)),MeansSim_mix(:,k-1),5,[0.8500 0.3250 0.0980],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(MeansSim_mix(:,k-1)) nanmean(MeansSim_mix(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(MeansSim_mix(:,k-1)) nanmedian(MeansSim_mix(:,k-1))],'k','linewidth',2)
end
title('G. Spatial random mixture'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Absolute agreement with original'); xlabel('Number of states')
subplot(2,4,8);
for k = 2:maxk
    swarmchart(k*ones(1,size(MeansSim_switch25,1)),MeansSim_gradient(:,k-1),5,[0.6350 0.0780 0.1840],'filled'); hold on;
    plot([k-0.25 k+0.25],[nanmean(MeansSim_gradient(:,k-1)) nanmean(MeansSim_gradient(:,k-1))],'k')
    plot([k-0.25 k+0.25],[nanmedian(MeansSim_gradient(:,k-1)) nanmedian(MeansSim_gradient(:,k-1))],'k','linewidth',2)
end
title('H. Spatial gradient'); axis([1.5 maxk+0.5 0 1]); set(gca,'xtick',2:maxk)
ylabel('Absolute agreement with original'); xlabel('Number of states')
print(gcf,'Results/HMM_Means_similarity_A1','-dpng','-r300');

figure
swarmchart(1*ones(1,size(MeansSim_switch50,1)),MeansSim_switch50(:,2),5,[0 0 1],'filled'); hold on;
plot([0.75 1.25],[nanmean(MeansSim_switch50(:,2)) nanmean(MeansSim_switch50(:,2))],'k'); plot([0.75 1.25],[nanmedian(MeansSim_switch50(:,2)) nanmedian(MeansSim_switch50(:,2))],'k','linewidth',2);
swarmchart(2*ones(1,size(MeansSim_switch25,1)),MeansSim_switch25(:,2),5,[0 1 1],'filled'); hold on;
plot([1.75 2.25],[nanmean(MeansSim_switch25(:,2)) nanmean(MeansSim_switch25(:,2))],'k'); plot([1.75 2.25],[nanmedian(MeansSim_switch25(:,2)) nanmedian(MeansSim_switch25(:,2))],'k','linewidth',2);
swarmchart(3*ones(1,size(MeansSim_switch10,1)),MeansSim_switch10(:,2),5,[0 0.4470 0.7410],'filled'); hold on;
plot([2.75 3.25],[nanmean(MeansSim_switch10(:,2)) nanmean(MeansSim_switch10(:,2))],'k'); plot([2.75 3.25],[nanmedian(MeansSim_switch10(:,2)) nanmedian(MeansSim_switch10(:,2))],'k','linewidth',2);
swarmchart(4*ones(1,size(MeansSim_max,1)),MeansSim_max(:,2),5,[0.4940 0.1840 0.5560],'filled'); hold on;
plot([3.75 4.25],[nanmean(MeansSim_max(:,2)) nanmean(MeansSim_max(:,2))],'k'); plot([3.75 4.25],[nanmedian(MeansSim_max(:,2)) nanmedian(MeansSim_max(:,2))],'k','linewidth',2);
swarmchart(5*ones(1,size(MeansSim_linear,1)),MeansSim_linear(:,2),5,[0 1 0],'filled'); hold on;
plot([4.75 5.25],[nanmean(MeansSim_linear(:,2)) nanmean(MeansSim_linear(:,2))],'k'); plot([4.75 5.25],[nanmedian(MeansSim_linear(:,2)) nanmedian(MeansSim_linear(:,2))],'k','linewidth',2);
swarmchart(6*ones(1,size(MeansSim_nonlinear,1)),MeansSim_nonlinear(:,2),5,[0.4660 0.6740 0.1880],'filled'); hold on;
plot([5.75 6.25],[nanmean(MeansSim_nonlinear(:,2)) nanmean(MeansSim_nonlinear(:,2))],'k'); plot([5.75 6.25],[nanmedian(MeansSim_nonlinear(:,2)) nanmedian(MeansSim_nonlinear(:,2))],'k','linewidth',2);
swarmchart(7*ones(1,size(MeansSim_mix,1)),MeansSim_mix(:,2),5,[0.8500 0.3250 0.0980],'filled'); hold on;
plot([6.75 7.25],[nanmean(MeansSim_mix(:,2)) nanmean(MeansSim_mix(:,2))],'k'); plot([6.75 7.25],[nanmedian(MeansSim_mix(:,2)) nanmedian(MeansSim_mix(:,2))],'k','linewidth',2);
swarmchart(8*ones(1,size(MeansSim_gradient,1)),MeansSim_gradient(:,2),5,[0.6350 0.0780 0.1840],'filled'); hold on;
plot([7.75 8.25],[nanmean(MeansSim_gradient(:,2)) nanmean(MeansSim_gradient(:,2))],'k'); plot([7.75 8.25],[nanmedian(MeansSim_gradient(:,2)) nanmedian(MeansSim_gradient(:,2))],'k','linewidth',2);
set(gca,'xtick',1:1:8,'xticklabel',{'Switch 50 TRs','Switch 25 TRs','Switch 10 TRs','Max switching','Linear additive coupling','Nonlinear multiplicative coupling','Spatial random mixture','Spatial gradient'},'FontSize',14)
ylabel('Absolute agreement with original','FontSize',14); title('HMM state mean similarity (based on 3 states)','FontSize',14)
print(gcf,'Results/HMM_Means_similarity_summary','-dpng','-r300');

figure
swarmchart(1*ones(1,size(TScorr,1)),TScorr(:,1),5,[0 0 1],'filled'); hold on;
plot([0.75 1.25],[nanmean(TScorr(:,1)) nanmean(TScorr(:,1))],'k'); plot([0.75 1.25],[nanmedian(TScorr(:,1)) nanmedian(TScorr(:,1))],'k','linewidth',2);
swarmchart(2*ones(1,size(TScorr,1)),TScorr(:,2),5,[0 1 1],'filled'); hold on;
plot([1.75 2.25],[nanmean(TScorr(:,2)) nanmean(TScorr(:,2))],'k'); plot([1.75 2.25],[nanmedian(TScorr(:,2)) nanmedian(TScorr(:,2))],'k','linewidth',2);
swarmchart(3*ones(1,size(TScorr,1)),TScorr(:,3),5,[0 0.4470 0.7410],'filled'); hold on;
plot([2.75 3.25],[nanmean(TScorr(:,3)) nanmean(TScorr(:,3))],'k'); plot([2.75 3.25],[nanmedian(TScorr(:,3)) nanmedian(TScorr(:,3))],'k','linewidth',2);
swarmchart(4*ones(1,size(TScorr,1)),TScorr(:,4),5,[0.4940 0.1840 0.5560],'filled'); hold on;
plot([3.75 4.25],[nanmean(TScorr(:,4)) nanmean(TScorr(:,4))],'k'); plot([3.75 4.25],[nanmedian(TScorr(:,4)) nanmedian(TScorr(:,4))],'k','linewidth',2);
swarmchart(5*ones(1,size(TScorr,1)),TScorr(:,5),5,[0 1 0],'filled'); hold on;
plot([4.75 5.25],[nanmean(TScorr(:,5)) nanmean(TScorr(:,5))],'k'); plot([4.75 5.25],[nanmedian(TScorr(:,5)) nanmedian(TScorr(:,5))],'k','linewidth',2);
swarmchart(6*ones(1,size(TScorr,1)),TScorr(:,6),5,[0.4660 0.6740 0.1880],'filled'); hold on;
plot([5.75 6.25],[nanmean(TScorr(:,6)) nanmean(TScorr(:,6))],'k'); plot([5.75 6.25],[nanmedian(TScorr(:,6)) nanmedian(TScorr(:,6))],'k','linewidth',2);
swarmchart(7*ones(1,size(TScorr,1)),TScorr(:,7),5,[0.8500 0.3250 0.0980],'filled'); hold on;
plot([6.75 7.25],[nanmean(TScorr(:,7)) nanmean(TScorr(:,7))],'k'); plot([6.75 7.25],[nanmedian(TScorr(:,7)) nanmedian(TScorr(:,7))],'k','linewidth',2);
swarmchart(8*ones(1,size(TScorr,1)),TScorr(:,8),5,[0.6350 0.0780 0.1840],'filled'); hold on;
plot([7.75 8.25],[nanmean(TScorr(:,8)) nanmean(TScorr(:,8))],'k'); plot([7.75 8.25],[nanmedian(TScorr(:,8)) nanmedian(TScorr(:,8))],'k','linewidth',2);
set(gca,'xtick',1:1:8,'xticklabel',{'Switch 50 TRs','Switch 25 TRs','Switch 10 TRs','Max switching','Linear additive coupling','Nonlinear multiplicative coupling','Spatial random mixture','Spatial gradient'},'FontSize',14)
ylabel('Correlation','FontSize',14); title('Timeseries similarity','FontSize',14)
print(gcf,'Results/TS_similarity','-dpng','-r300');



% Align means across participants
k = 3;
for s = setdiff(2:length(subs),[7 14])
    [~,assign] = matrixICC(Means_original_all{1,k-1},Means_original_all{s,k-1},1); Means_original_all(s,k-1) = {Means_original_all{s,k-1}(:,assign)};
    LifeTimes_original(s,:) = LifeTimes_original(s,assign); StateTransitions_original(s) = {StateTransitions_original{s}(assign,assign)};
    [~,assign] = matrixICC(Means_switch50_all{1,k-1},Means_switch50_all{s,k-1},1); Means_switch50_all(s,k-1) = {Means_switch50_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_switch25_all{1,k-1},Means_switch25_all{s,k-1},1); Means_switch25_all(s,k-1) = {Means_switch25_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_switch10_all{1,k-1},Means_switch10_all{s,k-1},1); Means_switch10_all(s,k-1) = {Means_switch10_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_max_all{1,k-1},Means_max_all{s,k-1},1); Means_max_all(s,k-1) = {Means_max_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_linear_all{1,k-1},Means_linear_all{s,k-1},1); Means_linear_all(s,k-1) = {Means_linear_all{s,k-1}(:,assign)};
    LifeTimes_linear(s,:) = LifeTimes_linear(s,assign); StateTransitions_linear(s) = {StateTransitions_linear{s}(assign,assign)};
    [~,assign] = matrixICC(Means_nonlinear_all{1,k-1},Means_nonlinear_all{s,k-1},1); Means_nonlinear_all(s,k-1) = {Means_nonlinear_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_mix_all{1,k-1},Means_mix_all{s,k-1},1); Means_mix_all(s,k-1) = {Means_mix_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_gradient_all{1,k-1},Means_gradient_all{s,k-1},1); Means_gradient_all(s,k-1) = {Means_gradient_all{s,k-1}(:,assign)};
end
M_original = []; M_switch50 = []; M_switch25 = []; M_switch10 = []; M_max = []; M_linear = []; M_nonlinear = []; M_mix = []; M_gradient = [];
for s = setdiff(1:length(subs),[7 14])
    M_original = cat(3,M_original,Means_original_all{s,k-1});
    M_switch50 = cat(3,M_switch50,Means_switch50_all{s,k-1});
    M_switch25 = cat(3,M_switch25,Means_switch25_all{s,k-1});
    M_switch10 = cat(3,M_switch10,Means_switch10_all{s,k-1});
    M_max = cat(3,M_max,Means_max_all{s,k-1});
    M_linear = cat(3,M_linear,Means_linear_all{s,k-1});
    M_nonlinear = cat(3,M_nonlinear,Means_nonlinear_all{s,k-1});
    M_mix = cat(3,M_mix,Means_mix_all{s,k-1});
    M_gradient = cat(3,M_gradient,Means_gradient_all{s,k-1});
end
figure
subplot(2,1,1);
swarmchart(1*ones(1,size(M_original,3)),squeeze(M_original(1,1,:)),5,'k','filled'); hold on;
plot([0.75 1.25],[nanmean(M_original(1,1,:),3) nanmean(M_original(1,1,:),3)],'k'); plot([0.75 1.25],[nanmedian(M_original(1,1,:),3) nanmedian(M_original(1,1,:),3)],'k','linewidth',2);
swarmchart(2*ones(1,size(M_original,3)),squeeze(M_original(2,1,:)),5,'k','filled'); hold on;
plot([1.75 2.25],[nanmean(M_original(2,1,:),3) nanmean(M_original(2,1,:),3)],'k'); plot([1.75 2.25],[nanmedian(M_original(2,1,:),3) nanmedian(M_original(2,1,:),3)],'k','linewidth',2);
swarmchart(3*ones(1,size(M_original,3)),squeeze(M_original(3,1,:)),5,'k','filled'); hold on;
plot([2.75 3.25],[nanmean(M_original(3,1,:),3) nanmean(M_original(3,1,:),3)],'k'); plot([2.75 3.25],[nanmedian(M_original(3,1,:),3) nanmedian(M_original(3,1,:),3)],'k','linewidth',2);
swarmchart(5*ones(1,size(M_original,3)),squeeze(M_original(1,2,:)),5,'k','filled'); hold on;
plot([4.75 5.25],[nanmean(M_original(1,2,:),3) nanmean(M_original(1,2,:),3)],'k'); plot([4.75 5.25],[nanmedian(M_original(1,2,:),3) nanmedian(M_original(1,2,:),3)],'k','linewidth',2);
swarmchart(6*ones(1,size(M_original,3)),squeeze(M_original(2,2,:)),5,'k','filled'); hold on;
plot([5.75 6.25],[nanmean(M_original(2,2,:),3) nanmean(M_original(2,2,:),3)],'k'); plot([5.75 6.25],[nanmedian(M_original(2,2,:),3) nanmedian(M_original(2,2,:),3)],'k','linewidth',2);
swarmchart(7*ones(1,size(M_original,3)),squeeze(M_original(3,2,:)),5,'k','filled'); hold on;
plot([6.75 7.25],[nanmean(M_original(3,2,:),3) nanmean(M_original(3,2,:),3)],'k'); plot([6.75 7.25],[nanmedian(M_original(3,2,:),3) nanmedian(M_original(3,2,:),3)],'k','linewidth',2);
swarmchart(9*ones(1,size(M_original,3)),squeeze(M_original(1,3,:)),5,'k','filled'); hold on;
plot([8.75 9.25],[nanmean(M_original(1,3,:),3) nanmean(M_original(1,3,:),3)],'k'); plot([8.75 9.25],[nanmedian(M_original(1,3,:),3) nanmedian(M_original(1,3,:),3)],'k','linewidth',2);
swarmchart(10*ones(1,size(M_original,3)),squeeze(M_original(2,3,:)),5,'k','filled'); hold on;
plot([9.75 10.25],[nanmean(M_original(2,3,:),3) nanmean(M_original(2,3,:),3)],'k'); plot([9.75 10.25],[nanmedian(M_original(2,3,:),3) nanmedian(M_original(2,3,:),3)],'k','linewidth',2);
swarmchart(11*ones(1,size(M_original,3)),squeeze(M_original(3,3,:)),5,'k','filled'); hold on;
plot([10.75 11.25],[nanmean(M_original(3,3,:),3) nanmean(M_original(3,3,:),3)],'k'); plot([10.75 11.25],[nanmedian(M_original(3,3,:),3) nanmedian(M_original(3,3,:),3)],'k','linewidth',2);
set(gca,'xtick',[1 2 3 5 6 7 9 10 11],'xticklabel',{'Network 1 - State 1','Network 2 - State 1','Overlap - State 1','Network 1 - State 2','Network 2 - State 2','Overlap - State 2','Network 1 - State 3','Network 2 - State 3','Overlap - State 3'})
ylabel('HMM Mean'); title('A. HMM means for original data'); hline(0,':k'); vline([4 8],'k-');
subplot(2,1,2);
swarmchart(1*ones(1,size(M_linear,3)),squeeze(M_linear(1,1,:)),5,'k','filled'); hold on;
plot([0.75 1.25],[nanmean(M_linear(1,1,:),3) nanmean(M_linear(1,1,:),3)],'k'); plot([0.75 1.25],[nanmedian(M_linear(1,1,:),3) nanmedian(M_linear(1,1,:),3)],'k','linewidth',2);
swarmchart(2*ones(1,size(M_linear,3)),squeeze(M_linear(2,1,:)),5,'k','filled'); hold on;
plot([1.75 2.25],[nanmean(M_linear(2,1,:),3) nanmean(M_linear(2,1,:),3)],'k'); plot([1.75 2.25],[nanmedian(M_linear(2,1,:),3) nanmedian(M_linear(2,1,:),3)],'k','linewidth',2);
swarmchart(3*ones(1,size(M_linear,3)),squeeze(M_linear(3,1,:)),5,'k','filled'); hold on;
plot([2.75 3.25],[nanmean(M_linear(3,1,:),3) nanmean(M_linear(3,1,:),3)],'k'); plot([2.75 3.25],[nanmedian(M_linear(3,1,:),3) nanmedian(M_linear(3,1,:),3)],'k','linewidth',2);
swarmchart(5*ones(1,size(M_linear,3)),squeeze(M_linear(1,2,:)),5,'k','filled'); hold on;
plot([4.75 5.25],[nanmean(M_linear(1,2,:),3) nanmean(M_linear(1,2,:),3)],'k'); plot([4.75 5.25],[nanmedian(M_linear(1,2,:),3) nanmedian(M_linear(1,2,:),3)],'k','linewidth',2);
swarmchart(6*ones(1,size(M_linear,3)),squeeze(M_linear(2,2,:)),5,'k','filled'); hold on;
plot([5.75 6.25],[nanmean(M_linear(2,2,:),3) nanmean(M_linear(2,2,:),3)],'k'); plot([5.75 6.25],[nanmedian(M_linear(2,2,:),3) nanmedian(M_linear(2,2,:),3)],'k','linewidth',2);
swarmchart(7*ones(1,size(M_linear,3)),squeeze(M_linear(3,2,:)),5,'k','filled'); hold on;
plot([6.75 7.25],[nanmean(M_linear(3,2,:),3) nanmean(M_linear(3,2,:),3)],'k'); plot([6.75 7.25],[nanmedian(M_linear(3,2,:),3) nanmedian(M_linear(3,2,:),3)],'k','linewidth',2);
swarmchart(9*ones(1,size(M_linear,3)),squeeze(M_linear(1,3,:)),5,'k','filled'); hold on;
plot([8.75 9.25],[nanmean(M_linear(1,3,:),3) nanmean(M_linear(1,3,:),3)],'k'); plot([8.75 9.25],[nanmedian(M_linear(1,3,:),3) nanmedian(M_linear(1,3,:),3)],'k','linewidth',2);
swarmchart(10*ones(1,size(M_linear,3)),squeeze(M_linear(2,3,:)),5,'k','filled'); hold on;
plot([9.75 10.25],[nanmean(M_linear(2,3,:),3) nanmean(M_linear(2,3,:),3)],'k'); plot([9.75 10.25],[nanmedian(M_linear(2,3,:),3) nanmedian(M_linear(2,3,:),3)],'k','linewidth',2);
swarmchart(11*ones(1,size(M_linear,3)),squeeze(M_linear(3,3,:)),5,'k','filled'); hold on;
plot([10.75 11.25],[nanmean(M_linear(3,3,:),3) nanmean(M_linear(3,3,:),3)],'k'); plot([10.75 11.25],[nanmedian(M_linear(3,3,:),3) nanmedian(M_linear(3,3,:),3)],'k','linewidth',2);
set(gca,'xtick',[1 2 3 5 6 7 9 10 11],'xticklabel',{'Network 1 - State 1','Network 2 - State 1','Overlap - State 1','Network 1 - State 2','Network 2 - State 2','Overlap - State 2','Network 1 - State 3','Network 2 - State 3','Overlap - State 3'})
ylabel('HMM Mean'); title('B. HMM means for simulated linear additive coupling data'); hline(0,':k'); vline([4 8],'k-');
print(gcf,'Results/HMM_mean_original_linear_k3','-dpng','-r300');

% Align means across participants
k = 2;
for s = setdiff(2:length(subs),[7 14])
    [~,assign] = matrixICC(Means_original_all{1,k-1},Means_original_all{s,k-1},1); Means_original_all(s,k-1) = {Means_original_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_switch50_all{1,k-1},Means_switch50_all{s,k-1},1); Means_switch50_all(s,k-1) = {Means_switch50_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_switch25_all{1,k-1},Means_switch25_all{s,k-1},1); Means_switch25_all(s,k-1) = {Means_switch25_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_switch10_all{1,k-1},Means_switch10_all{s,k-1},1); Means_switch10_all(s,k-1) = {Means_switch10_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_max_all{1,k-1},Means_max_all{s,k-1},1); Means_max_all(s,k-1) = {Means_max_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_linear_all{1,k-1},Means_linear_all{s,k-1},1); Means_linear_all(s,k-1) = {Means_linear_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_nonlinear_all{1,k-1},Means_nonlinear_all{s,k-1},1); Means_nonlinear_all(s,k-1) = {Means_nonlinear_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_mix_all{1,k-1},Means_mix_all{s,k-1},1); Means_mix_all(s,k-1) = {Means_mix_all{s,k-1}(:,assign)};
    [~,assign] = matrixICC(Means_gradient_all{1,k-1},Means_gradient_all{s,k-1},1); Means_gradient_all(s,k-1) = {Means_gradient_all{s,k-1}(:,assign)};
end
M_original = []; M_switch50 = []; M_switch25 = []; M_switch10 = []; M_max = []; M_linear = []; M_nonlinear = []; M_mix = []; M_gradient = [];
for s = setdiff(1:length(subs),[7 14])
    M_original = cat(3,M_original,Means_original_all{s,k-1});
    M_switch50 = cat(3,M_switch50,Means_switch50_all{s,k-1});
    M_switch25 = cat(3,M_switch25,Means_switch25_all{s,k-1});
    M_switch10 = cat(3,M_switch10,Means_switch10_all{s,k-1});
    M_max = cat(3,M_max,Means_max_all{s,k-1});
    M_linear = cat(3,M_linear,Means_linear_all{s,k-1});
    M_nonlinear = cat(3,M_nonlinear,Means_nonlinear_all{s,k-1});
    M_mix = cat(3,M_mix,Means_mix_all{s,k-1});
    M_gradient = cat(3,M_gradient,Means_gradient_all{s,k-1});
end
figure
subplot(2,1,1);
swarmchart(1*ones(1,size(M_original,3)),squeeze(M_original(1,1,:)),5,'k','filled'); hold on;
plot([0.75 1.25],[nanmean(M_original(1,1,:),3) nanmean(M_original(1,1,:),3)],'k'); plot([0.75 1.25],[nanmedian(M_original(1,1,:),3) nanmedian(M_original(1,1,:),3)],'k','linewidth',2);
swarmchart(2*ones(1,size(M_original,3)),squeeze(M_original(2,1,:)),5,'k','filled'); hold on;
plot([1.75 2.25],[nanmean(M_original(2,1,:),3) nanmean(M_original(2,1,:),3)],'k'); plot([1.75 2.25],[nanmedian(M_original(2,1,:),3) nanmedian(M_original(2,1,:),3)],'k','linewidth',2);
swarmchart(3*ones(1,size(M_original,3)),squeeze(M_original(3,1,:)),5,'k','filled'); hold on;
plot([2.75 3.25],[nanmean(M_original(3,1,:),3) nanmean(M_original(3,1,:),3)],'k'); plot([2.75 3.25],[nanmedian(M_original(3,1,:),3) nanmedian(M_original(3,1,:),3)],'k','linewidth',2);
swarmchart(5*ones(1,size(M_original,3)),squeeze(M_original(1,2,:)),5,'k','filled'); hold on;
plot([4.75 5.25],[nanmean(M_original(1,2,:),3) nanmean(M_original(1,2,:),3)],'k'); plot([4.75 5.25],[nanmedian(M_original(1,2,:),3) nanmedian(M_original(1,2,:),3)],'k','linewidth',2);
swarmchart(6*ones(1,size(M_original,3)),squeeze(M_original(2,2,:)),5,'k','filled'); hold on;
plot([5.75 6.25],[nanmean(M_original(2,2,:),3) nanmean(M_original(2,2,:),3)],'k'); plot([5.75 6.25],[nanmedian(M_original(2,2,:),3) nanmedian(M_original(2,2,:),3)],'k','linewidth',2);
swarmchart(7*ones(1,size(M_original,3)),squeeze(M_original(3,2,:)),5,'k','filled'); hold on;
plot([6.75 7.25],[nanmean(M_original(3,2,:),3) nanmean(M_original(3,2,:),3)],'k'); plot([6.75 7.25],[nanmedian(M_original(3,2,:),3) nanmedian(M_original(3,2,:),3)],'k','linewidth',2);
set(gca,'xtick',[1 2 3 5 6 7],'xticklabel',{'Network 1 - State 1','Network 2 - State 1','Overlap - State 1','Network 1 - State 2','Network 2 - State 2','Overlap - State 2'})
ylabel('HMM Mean'); title('A. HMM means for original data'); hline(0,':k'); vline([4],'k-');
subplot(2,1,2);
swarmchart(1*ones(1,size(M_linear,3)),squeeze(M_linear(1,1,:)),5,'k','filled'); hold on;
plot([0.75 1.25],[nanmean(M_linear(1,1,:),3) nanmean(M_linear(1,1,:),3)],'k'); plot([0.75 1.25],[nanmedian(M_linear(1,1,:),3) nanmedian(M_linear(1,1,:),3)],'k','linewidth',2);
swarmchart(2*ones(1,size(M_linear,3)),squeeze(M_linear(2,1,:)),5,'k','filled'); hold on;
plot([1.75 2.25],[nanmean(M_linear(2,1,:),3) nanmean(M_linear(2,1,:),3)],'k'); plot([1.75 2.25],[nanmedian(M_linear(2,1,:),3) nanmedian(M_linear(2,1,:),3)],'k','linewidth',2);
swarmchart(3*ones(1,size(M_linear,3)),squeeze(M_linear(3,1,:)),5,'k','filled'); hold on;
plot([2.75 3.25],[nanmean(M_linear(3,1,:),3) nanmean(M_linear(3,1,:),3)],'k'); plot([2.75 3.25],[nanmedian(M_linear(3,1,:),3) nanmedian(M_linear(3,1,:),3)],'k','linewidth',2);
swarmchart(5*ones(1,size(M_linear,3)),squeeze(M_linear(1,2,:)),5,'k','filled'); hold on;
plot([4.75 5.25],[nanmean(M_linear(1,2,:),3) nanmean(M_linear(1,2,:),3)],'k'); plot([4.75 5.25],[nanmedian(M_linear(1,2,:),3) nanmedian(M_linear(1,2,:),3)],'k','linewidth',2);
swarmchart(6*ones(1,size(M_linear,3)),squeeze(M_linear(2,2,:)),5,'k','filled'); hold on;
plot([5.75 6.25],[nanmean(M_linear(2,2,:),3) nanmean(M_linear(2,2,:),3)],'k'); plot([5.75 6.25],[nanmedian(M_linear(2,2,:),3) nanmedian(M_linear(2,2,:),3)],'k','linewidth',2);
swarmchart(7*ones(1,size(M_linear,3)),squeeze(M_linear(3,2,:)),5,'k','filled'); hold on;
plot([6.75 7.25],[nanmean(M_linear(3,2,:),3) nanmean(M_linear(3,2,:),3)],'k'); plot([6.75 7.25],[nanmedian(M_linear(3,2,:),3) nanmedian(M_linear(3,2,:),3)],'k','linewidth',2);
set(gca,'xtick',[1 2 3 5 6 7],'xticklabel',{'Network 1 - State 1','Network 2 - State 1','Overlap - State 1','Network 1 - State 2','Network 2 - State 2','Overlap - State 2'})
ylabel('HMM Mean'); title('B. HMM means for simulated linear additive coupling data'); hline(0,':k'); vline([4],'k-');
print(gcf,'Results/HMM_mean_original_linear_k2','-dpng','-r300');

% Statistics
nanmedian(MeansSim_linear)
nanmean(MeansSim_linear)
[H,P,CI,STATS] = ttest2(MeansSim_linear(:,2),MeansSim_switch10(:,2))
[H,P,CI,STATS] = ttest2(MeansSim_linear(:,2),MeansSim_linear(:,1))

% Save output
clear r rep1 rep2 reps s Scans tp v window D d assign i I I7 I7new I8 I8new k map7 map8
save('Results/HMM_output.mat');

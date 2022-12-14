clear all; close all; clc

addpath(genpath('/scratch/janine.bijsterbosch/HCP/Matlab/'),'-END')

subs = [105923 114823 125525 130518 137128 144226 146129 158035 169343 177746 185442 192439 195041 200614 204521 562345 601127 783462 859671 861456];
scans = {'3T/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii',...
    '3T/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii',...
    '3T/rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii',...
    '3T/rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii',...
    'retest/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii',...
    'retest/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii',...
    'retest/rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii',...
    'retest/rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii',...
    '7T/rfMRI_REST1_7T_Atlas_MSMAll_hp2000_clean.dtseries.nii',...
    '7T/rfMRI_REST2_7T_Atlas_MSMAll_hp2000_clean.dtseries.nii',...
    '7T/rfMRI_REST3_7T_Atlas_MSMAll_hp2000_clean.dtseries.nii',...
    '7T/rfMRI_REST4_7T_Atlas_MSMAll_hp2000_clean.dtseries.nii'};

load PFM_keep.mat arenans;

%% [steve] janine's follow-up on tendency for overlap regions to have negative skew
subject_maps_separate_runs = ft_read_cifti('Results/Maps/subject_maps_separate_runs.dtseries.nii');
subject_maps_separate_runs = subject_maps_separate_runs.dtseries;
subject_maps_separate_runs(arenans,:) = [];
SumMapAll = zeros(91282,length(subs));
SkewMapAll = zeros(91282,length(subs));
for s = 1:length(subs)
    fprintf('Running subject %d\n',s)
    M = subject_maps_separate_runs(:,(s-1)*12+1:s*12);
    SumMapAll(:,s) = sum(M,2); clear M
    SK = zeros(91282,length(scans));
    for scan = 1:length(scans)
        D = ft_read_cifti(sprintf('/scratch/janine.bijsterbosch/HCP/%d/%s',subs(s),scans{scan}));
        D = D.dtseries;
        D(isnan(D(:,1))==1,:) = [];
        SK(:,scan) = skewness(D,1,2); clear D
    end
    SkewMapAll(:,s) = mean(SK,2); clear SK
end
save('Skew_Sum_test.mat','SkewMapAll','SumMapAll')

%% Run locally after downloading Skew_Sum_test.mat from CHPC


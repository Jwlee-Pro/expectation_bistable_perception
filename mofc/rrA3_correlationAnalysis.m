clear all; close all; clc ;


ValSub =  1:21; 
code_loc = '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception'; 

for iSub = ValSub
    load([code_loc '/data/wholebrain/corr_sub' num2str(iSub) '.mat'],'Corr_simple', 'Corr_part')
    cMat_simple(:,iSub) = Corr_simple.sfm_corr; 
    cMat_part(:,iSub) = Corr_part.sfm_corr; 
end

crits = ones(length(GrayCoord),1); 
nifti_mapping(mean(cMat_simple.^2,2), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_simple')
nifti_mapping(mean(cMat_part.^2,2), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_part')



%% Compare DV, Choice, and expectation 
for iSub = ValSub
    iSub
    load([code_loc '/data/wholebrain/corr_sub' num2str(iSub) '.mat'],'Corr_simple', 'Corr_part')
    cMat_simple(:,iSub) = Corr_simple.sfm_corr; 
    cMat_part(:,iSub) = Corr_part.sfm_corr; 
    
    load([code_loc '/data/wholebrain/DV/corr_sub' num2str(iSub) '.mat'],'Corr_simple')
    dv_simple(:,iSub) = Corr_simple.sfm_corr; 
    
    load([code_loc '/data/wholebrain/Decision/corr_sub' num2str(iSub) '.mat'],'Corr_simple')
    dec_simple(:,iSub) = Corr_simple.sfm_corr; 
end

[h1,p,ci,stats1] = ttest((cMat_simple.^2 - dv_simple.^2)');
nifti_mapping(mean((cMat_simple.^2 - dv_simple.^2)'), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'dv_compare')
[h2,p,ci,stats2] = ttest((cMat_simple.^2 - dec_simple.^2)');
nifti_mapping(mean((cMat_simple.^2 - dec_simple.^2)'), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'dec_compare')


% Apply multiple criterions 
temp = mean(cMat_part.^2,2); 
nifti_mapping(temp, temp>quantile(temp,0.7) & h1'==1 & stats1.tstat'>0  & h2'==1 & stats2.tstat'>0, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_part_multiplecrits')




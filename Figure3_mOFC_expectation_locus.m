%% Locus of expectation in medial OFC
% Coded by Joonwon Lee

clear all; close all; clc; 


%% [0] Module details
% 1. Whole brain GLM
% 2. ROI selection, drawing clusters 
% 3. SVR train-test (multivariate pattern analysis) 
% 4. Validation in trial-to-trial level 
% 5. Validation lawful relationship 
runModules = {};
runModules.compute_corr = 1; 
runModules.criterion_filter = 1; 


%% [1] Preparation 
global gParam time_fmri
gParam.TR = 2.4;                % TR
gParam.HD = 2;                  % Hemodynamic delay (for HD shift and population decoding) 
gParam.EventWindow = 34;        % +-20sec window for Transition event 
gParam.FigOn = 1;
gParam.TrCrit = [2 250 2];     % Transition trial selection criterion (pre-stable, unstable, post-stable duration)
time_fmri = (1:125)*gParam.TR;
gParam.time_fmri_rs = time_fmri(1):0.1:time_fmri(end); 
cmap_choice=[102 204 255;204 204 204; 255 204 102]./255;
Common_fontsize = 10; 
ValSub = 1:21; 


code_loc = '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception'; 


% Data loading information & toolbox setting
fprintf('Data loading information & toolbox setting.... \n'); 
addpath(genpath([code_loc '/src/packages/']));
addpath([code_loc '/src/library/']); 

% Load behavioral data (model regressors)
fprintf('Load behavioral data.... \n'); 
ParamSetting_fmri;
load(['/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/' ...
    'Paper Writing/LeeEtal_ManuSfM/2017_M7/MTanalysis_deep/' ...
    'LeeLee_PriorEffect_workspace_PredCo20171127.mat'],'matSFM','matMimic','matSFM_PredCo'); 
load([code_loc '/data/neural/mDataMat_all']); 

Percept_all = cell(1,21); Percept_all_mimic = cell(1,21); 
PE_all = cell(1,21); PE_all_mimic = cell(1,21); 
U_all = cell(1,21); H_all = cell(1,21);
DV_all = cell(1,21); DV_all_mimic = cell(1,21);
nSess = [ones(1,14) ones(1,7)*4];
for iSub = ValSub
    for iS = 1:4
        Percept_all{iSub} = [Percept_all{iSub} mDataMat.Decision_conv{iSub,iS}]; 
        DV_all{iSub} = [DV_all{iSub} mDataMat.DV_conv{iSub,iS}]; 
        PE_all{iSub} = [PE_all{iSub} mDataMat.PE_conv{iSub,iS}];         
        H_all{iSub} = [H_all{iSub} mDataMat.H_conv{iSub,iS}];         
        U_all{iSub} = [U_all{iSub} mDataMat.U_conv{iSub,iS}];         
    end
    
    for iS = 1:nSess(iSub)
        Percept_all_mimic{iSub} = [Percept_all_mimic{iSub} matMimic{iS,iSub}.result_rs.Percept]; 
        PE_all_mimic{iSub} = [PE_all_mimic{iSub} matMimic{iS,iSub}.result_rs.matPrior];         
    end
end

% Load gray matter coordinates
ROI_Gray = load_nii([code_loc '/data/neural/grey_dimMatched.nii']);
[nx,ny,nz,nt] = size(ROI_Gray.img);
GrayCoord = []; 
for ix = 1:nx
    for iy = 1:ny
        for iz = 1:nz
            if ROI_Gray.img(ix,iy,iz) > 0.01
                GrayCoord = [GrayCoord; [ix iy iz]];
            end
        end
    end
end



%% Compute correlations 
if runModules.compute_corr == 1 
    for iSub = ValSub
        fprintf(['Computing sub-' num2str(iSub) '\n']); 
        load([code_loc '/data/wholebrain/wb_bold_sub' num2str(iSub) '.mat'],'Nifti_AlliS_filt', 'Nifti_AlliS_rep_filt')
        % temp_pe = diff(PE_all{iSub}); 
        % temp_c = diff(Percept_all{iSub}); 
        % temp_dv = diff(DV_all{iSub}); 
        % for iVox = 1:length(GrayCoord)
        %     nanInd = ~isnan(Nifti_AlliS_filt(iVox,:) + PE_all{iSub} + Percept_all{iSub} + DV_all{iSub}); 
        %     nanInd_diff = ~isnan(diff(Nifti_AlliS_filt(iVox,:)) + diff(PE_all{iSub}) + diff(Percept_all{iSub}) + diff(DV_all{iSub})); 
        %     if sum(nanInd) > 5
        %         [corrMat.pe(iVox), pvalMat.pe(iVox)] = corr(Nifti_AlliS_filt(iVox,nanInd)', PE_all{iSub}(nanInd)'); 
        %         [corrMat.dv(iVox), pvalMat.dv(iVox)] = corr(Nifti_AlliS_filt(iVox,nanInd)', DV_all{iSub}(nanInd)'); 
        %         [corrMat.c(iVox), pvalMat.c(iVox)] = corr(Nifti_AlliS_filt(iVox,nanInd)', Percept_all{iSub}(nanInd)'); 
        % 
        %         [corrMat.pe_c(iVox), temp] = partcorr(Nifti_AlliS_filt(iVox,nanInd)', PE_all{iSub}(nanInd)', Percept_all{iSub}(nanInd)'); 
        %         pvalMat.pe_c(iVox) = temp.p; 
        %         [corrMat.c_pe(iVox), temp] = partcorr(Nifti_AlliS_filt(iVox,nanInd)', Percept_all{iSub}(nanInd)', PE_all{iSub}(nanInd)'); 
        %         pvalMat.c_pe(iVox) = temp.p; 
        % 
        %         [corrMat.pe_dv(iVox), temp] = partcorr(Nifti_AlliS_filt(iVox,nanInd)', PE_all{iSub}(nanInd)', DV_all{iSub}(nanInd)'); 
        %         pvalMat.pe_dv(iVox) = temp.p; 
        %         [corrMat.dv_pe(iVox), temp] = partcorr(Nifti_AlliS_filt(iVox,nanInd)', DV_all{iSub}(nanInd)', PE_all{iSub}(nanInd)'); 
        %         pvalMat.dv_pe(iVox) = temp.p; 
        % 
        %         % Diff
        %         tempv = diff(Nifti_AlliS_filt(iVox,:)); 
        %         [corrMat.diff_pe(iVox), pvalMat.diff_pe(iVox)] = corr(tempv(nanInd_diff)', temp_pe(nanInd_diff)'); 
        %         [corrMat.diff_dv(iVox), pvalMat.diff_dv(iVox)] = corr(tempv(nanInd_diff)', temp_dv(nanInd_diff)'); 
        %         [corrMat.diff_c(iVox), pvalMat.diff_c(iVox)] = corr(tempv(nanInd_diff)', temp_c(nanInd_diff)'); 
        %     else
        %         corrMat.pe(iVox) = nan; pvalMat.pe(iVox) = nan; 
        %         corrMat.dv(iVox) = nan; pvalMat.dv(iVox) = nan; 
        %         corrMat.c(iVox) = nan; pvalMat.c(iVox) = nan; 
        % 
        %         corrMat.pe_c(iVox) = nan; pvalMat.pe_c(iVox) = nan; 
        %         corrMat.c_pe(iVox) = nan; pvalMat.c_pe(iVox) = nan; 
        %         corrMat.pe_dv(iVox) = nan; pvalMat.pe_dv(iVox) = nan; 
        %         corrMat.dv_pe(iVox) = nan; pvalMat.dv_pe(iVox) = nan; 
        % 
        %         corrMat.diff_pe(iVox) = nan; pvalMat.diff_pe(iVox) = nan; 
        %         corrMat.diff_dv(iVox) = nan; pvalMat.diff_dv(iVox) = nan; 
        %         corrMat.diff_c(iVox) = nan; pvalMat.diff_c(iVox) = nan; 
        % 
        %     end
        % end
        % save([code_loc '/data/wholebrain/corr_sub' num2str(iSub) '.mat'],'corrMat', 'pvalMat')


        % Replay 
        temp_pe = diff(PE_all_mimic{iSub}); 
        temp_c = diff(Percept_all_mimic{iSub}); 
        for iVox = 1:length(GrayCoord)
            nanInd = ~isnan(Nifti_AlliS_rep_filt(iVox,:) + PE_all_mimic{iSub} + Percept_all_mimic{iSub}); 
            nanInd_diff = ~isnan(diff(Nifti_AlliS_rep_filt(iVox,:)) + diff(PE_all_mimic{iSub}) + diff(Percept_all_mimic{iSub})); 
            if sum(nanInd) > 5
                [corrMat.pe(iVox), pvalMat.pe(iVox)] = corr(Nifti_AlliS_rep_filt(iVox,nanInd)', PE_all_mimic{iSub}(nanInd)'); 
                % [corrMat.dv(iVox), pvalMat.dv(iVox)] = corr(Nifti_AlliS_rep_filt(iVox,nanInd)', DV_all{iSub}(nanInd)'); 
                [corrMat.c(iVox), pvalMat.c(iVox)] = corr(Nifti_AlliS_rep_filt(iVox,nanInd)', Percept_all_mimic{iSub}(nanInd)'); 

                [corrMat.pe_c(iVox), temp] = partcorr(Nifti_AlliS_rep_filt(iVox,nanInd)', PE_all_mimic{iSub}(nanInd)', Percept_all_mimic{iSub}(nanInd)'); 
                pvalMat.pe_c(iVox) = temp.p; 
                [corrMat.c_pe(iVox), temp] = partcorr(Nifti_AlliS_rep_filt(iVox,nanInd)', Percept_all_mimic{iSub}(nanInd)', PE_all_mimic{iSub}(nanInd)'); 
                pvalMat.c_pe(iVox) = temp.p; 

                % [corrMat.pe_dv(iVox), temp] = partcorr(Nifti_AlliS_rep_filt(iVox,nanInd)', PE_all_mimic{iSub}(nanInd)', DV_all{iSub}(nanInd)'); 
                % pvalMat.pe_dv(iVox) = temp.p; 
                % [corrMat.dv_pe(iVox), temp] = partcorr(Nifti_AlliS_rep_filt(iVox,nanInd)', DV_all{iSub}(nanInd)', PE_all_mimic{iSub}(nanInd)'); 
                % pvalMat.dv_pe(iVox) = temp.p; 
                
                % Diff
                tempv = diff(Nifti_AlliS_rep_filt(iVox,:)); 
                [corrMat.diff_pe(iVox), pvalMat.diff_pe(iVox)] = corr(tempv(nanInd_diff)', temp_pe(nanInd_diff)'); 
                % [corrMat.diff_dv(iVox), pvalMat.diff_dv(iVox)] = corr(tempv(nanInd_diff)', temp_dv(nanInd_diff)'); 
                [corrMat.diff_c(iVox), pvalMat.diff_c(iVox)] = corr(tempv(nanInd_diff)', temp_c(nanInd_diff)'); 
            else
                corrMat.pe(iVox) = nan; pvalMat.pe(iVox) = nan; 
                % corrMat.dv(iVox) = nan; pvalMat.dv(iVox) = nan; 
                corrMat.c(iVox) = nan; pvalMat.c(iVox) = nan; 
                
                corrMat.pe_c(iVox) = nan; pvalMat.pe_c(iVox) = nan; 
                corrMat.c_pe(iVox) = nan; pvalMat.c_pe(iVox) = nan; 
                % corrMat.pe_dv(iVox) = nan; pvalMat.pe_dv(iVox) = nan; 
                % corrMat.dv_pe(iVox) = nan; pvalMat.dv_pe(iVox) = nan; 
                
                corrMat.diff_pe(iVox) = nan; pvalMat.diff_pe(iVox) = nan; 
                % corrMat.diff_dv(iVox) = nan; pvalMat.diff_dv(iVox) = nan; 
                corrMat.diff_c(iVox) = nan; pvalMat.diff_c(iVox) = nan; 
                
            end
        end
        save([code_loc '/data/wholebrain/corr_replay_sub' num2str(iSub) '.mat'],'corrMat', 'pvalMat')
    end
end
asdf

%% Find voxels satisfying criterions
if runModules.criterion_filter == 1 
    crits = ones(length(GrayCoord),1); 
    mat_pe = []; 
    mat_dv = []; 
    mat_c = []; 
    pmat_pe_c = []; 
    pmat_pe_dv = []; 
    pmat_c_pe = []; 
    diff_pe = []; 
    
    for iSub = ValSub
        load([code_loc '/data/wholebrain/corr_sub' num2str(iSub) '.mat'],'corrMat', 'pvalMat')
        mat_pe = [mat_pe; corrMat.pe];
        mat_dv = [mat_dv; corrMat.dv];
        mat_c = [mat_c; corrMat.c];
        
        pmat_pe_c = [pmat_pe_c; corrMat.pe_c];
        pmat_pe_dv = [pmat_pe_dv; corrMat.pe_dv];
        pmat_c_pe = [pmat_c_pe; corrMat.c_pe];
        diff_pe = [diff_pe; corrMat.diff_pe]; 
    end
    
    % Nifit-mapping & save files 
    nifti_mapping(mean(mat_pe.^2,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe')
    nifti_mapping(mean(mat_dv.^2,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_dv')
    nifti_mapping(mean(mat_c.^2,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_c')
    
%     nifti_mapping(mean(mat_pe_c.^2,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe_c')
%     nifti_mapping(mean(mat_pe_dv.^2,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe_dv')
%     nifti_mapping(mean(mat_c_pe.^2,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_c_pe')
    
    for iSub = ValSub
        nifti_mapping(mat_pe(iSub,:).^2, crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural/individual', ['varmap_pe_sub' num2str(iSub)])
    end
    
    % Relative information
    nifti_mapping(mean(mat_pe.^2,1), (mean(mat_pe.^2,1) - mean(mat_c.^2,1))>0, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe-c')
    nifti_mapping(mean(mat_pe.^2,1), (mean(mat_pe.^2,1) - mean(mat_dv.^2,1))>0, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe-dv')
    nifti_mapping(mean(mat_pe.^2,1), (mean(mat_pe.^2,1) - mean(mat_dv.^2,1))>0 & (mean(mat_pe.^2,1) - mean(mat_c.^2,1))>0, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe_crit1')

    nifti_mapping(mean(mat_pe.^2,1), (mean(mat_pe.^2,1) - mean(mat_c.^2,1))>0 & (mean(mat_pe.^2,1) > quantile(mean(mat_pe.^2,1),0.95)), ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe_crit1_cut1')
    nifti_mapping(mean(mat_pe.^2,1), (mean(mat_pe.^2,1) - mean(mat_dv.^2,1))>0 & (mean(mat_pe.^2,1) - mean(mat_c.^2,1))>0 & (mean(mat_pe.^2,1) > quantile(mean(mat_pe.^2,1),0.95)), ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe_crit1_cut2')
    
end



%% Step #1
% correlated to expectation
nifti_mapping(nanmean(mat_pe.^2,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe')

figure(100); clf; 
SP = subplot(1,1,1); cla; hold on; 
hist(nanmean(mat_pe.^2,1),40)
[N,X] = hist(nanmean(mat_pe.^2,1),40); 
Bh = bar(X,N,'facecolor',[0 0 0]+0.7);

vv = quantile(mean(mat_pe.^2,1),0.95); 
nifti_mapping(mean(mat_pe.^2,1), mean(mat_pe.^2,1)>vv, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe_crit')
nifti_mapping(mean(mat_pe.^2,1), (mean(mat_pe.^2,1) - mean(mat_c.^2,1))>0 & mean(mat_pe.^2,1)>vv, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe-c_crit')

nifti_mapping(mean(mat_pe.^2,1), mean(mat_pe.^2,1)>vv, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe_crit')
nifti_mapping(mean(mat_pe.^2,1), (mean(mat_pe.^2,1) - mean(mat_dv.^2,1))>0 & mean(mat_pe.^2,1)>vv, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe-dv_crit')

nifti_mapping(mean(mat_pe.^2,1), mean(mat_pe.^2,1)>vv, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe_crit')
nifti_mapping(mean(mat_pe.^2,1), (mean(mat_pe.^2,1) - mean(mat_dv.^2,1))>0 & mean(mat_pe.^2,1)>vv, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'pmat_pe_c_crit')


nifti_mapping(mean(diff_pe,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'diffvar_pe_crit')

% Invariant to current choice

[x1, x2] = ttest(pmat_c_pe.*sign(mat_pe));

nifti_mapping(mean(mat_pe.^2,1), (mean(mat_pe.^2,1) - mean(mat_c.^2,1))>0 & mean(mat_pe.^2,1)>vv, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe-c_crit')

% Sign-matched part correlation (c - pe)
nifti_mapping(mean(pmat_c_pe.*sign(mat_pe),1), (mean(mat_pe.^2,1) - mean(mat_c.^2,1))>0 & mean(mat_pe.^2,1)>vv & x1==0, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe-c_cv_crit')



%% Step #2 Replay
if runModules.criterion_filter == 1 
    crits = ones(length(GrayCoord),1); 
    mat_pe = []; 
    mat_c = []; 
    pmat_pe_c = []; 
    pmat_c_pe = []; 
    diff_pe = []; 
    
    for iSub = ValSub
        load([code_loc '/data/wholebrain/corr_replay_sub' num2str(iSub) '.mat'],'corrMat', 'pvalMat')
        mat_pe = [mat_pe; corrMat.pe];
        mat_c = [mat_c; corrMat.c];
        
        pmat_pe_c = [pmat_pe_c; corrMat.pe_c];
        pmat_c_pe = [pmat_c_pe; corrMat.c_pe];
        diff_pe = [diff_pe; corrMat.diff_pe]; 
    end
    
    % Nifit-mapping & save files 
    nifti_mapping(nanmean(mat_pe.^2,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_rep_pe')
    nifti_mapping(mean(mat_c.^2,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_rep_c')
    
    % Relative information
    nifti_mapping(mean(mat_pe.^2,1), (mean(mat_pe.^2,1) - mean(mat_c.^2,1))>0, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_rep_pe-c')

    vv = quantile(mean(mat_pe.^2,1),0.95); 
    nifti_mapping(mean(mat_pe.^2,1), mean(mat_pe.^2,1)>vv, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe_crit')
    nifti_mapping(mean(mat_pe.^2,1), (mean(mat_pe.^2,1) - mean(mat_c.^2,1))>0 & mean(mat_pe.^2,1)>vv, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_rep_pe-c_crit')
end

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

        for iVox = 1:length(GrayCoord)
            nanInd = ~isnan(Nifti_AlliS_filt(iVox,:) + PE_all{iSub} + Percept_all{iSub} + DV_all{iSub}); 
            if sum(nanInd) > 5
                [corrMat.pe(iVox), pvalMat.pe(iVox)] = corr(Nifti_AlliS_filt(iVox,nanInd)', PE_all{iSub}(nanInd)'); 
                [corrMat.dv(iVox), pvalMat.dv(iVox)] = corr(Nifti_AlliS_filt(iVox,nanInd)', DV_all{iSub}(nanInd)'); 
                [corrMat.c(iVox), pvalMat.c(iVox)] = corr(Nifti_AlliS_filt(iVox,nanInd)', Percept_all{iSub}(nanInd)'); 

                [corrMat.pe_c(iVox), temp] = partcorr(Nifti_AlliS_filt(iVox,nanInd)', PE_all{iSub}(nanInd)', Percept_all{iSub}(nanInd)'); 
                pvalMat.pe_c(iVox) = temp.p; 
                [corrMat.c_pe(iVox), temp] = partcorr(Nifti_AlliS_filt(iVox,nanInd)', Percept_all{iSub}(nanInd)', PE_all{iSub}(nanInd)'); 
                pvalMat.c_pe(iVox) = temp.p; 

                [corrMat.pe_dv(iVox), temp] = partcorr(Nifti_AlliS_filt(iVox,nanInd)', PE_all{iSub}(nanInd)', DV_all{iSub}(nanInd)'); 
                pvalMat.pe_dv(iVox) = temp.p; 
                [corrMat.dv_pe(iVox), temp] = partcorr(Nifti_AlliS_filt(iVox,nanInd)', DV_all{iSub}(nanInd)', PE_all{iSub}(nanInd)'); 
                pvalMat.dv_pe(iVox) = temp.p; 
            end
        end
        save([code_loc '/data/wholebrain/corr_sub' num2str(iSub) '.mat'],'corrMat', 'pvalMat')
    end
end

%% Find voxels satisfying criterions
if runModules.criterion_filter == 1 
    crits = ones(length(GrayCoord),1); 
    mat_pe = []; 
    mat_dv = []; 
    mat_c = []; 
    for iSub = ValSub
        load([code_loc '/data/wholebrain/corr_sub' num2str(iSub) '.mat'],'corrMat', 'pvalMat')
        mat_pe = [mat_pe; corrMat.pe];
        mat_dv = [mat_dv; corrMat.dv];
        mat_c = [mat_c; corrMat.c];
    end
    
    % Nifit-mapping & save files 
    nifti_mapping(mean(mat_pe.^2,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_pe')
    nifti_mapping(mean(mat_dv.^2,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_dv')
    nifti_mapping(mean(mat_c.^2,1), crits, ROI_Gray, GrayCoord, '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception/data/neural', 'varmap_c')
end



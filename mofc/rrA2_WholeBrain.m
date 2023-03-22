%% Locus of expectation in medial OFC
% Coded by Joonwon Lee

clear all; close all; clc; 


%% [0] Module details
% 1. Whole brain GLM
% 2. ROI selection, drawing clusters 
% 3. SVR train-test (multivariate pattern analysis) 
% 4. Validation in trial-to-trial level 
% 5. Validation lawful relationship 



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


code_loc = '/Volumes/ROOT/CSNL_temp/JWL/expectation_bistable_perception'; 


% Data loading information & toolbox setting
fprintf('Data loading information & toolbox setting.... \n'); 
addpath(genpath([code_loc '/src/packages/']));
addpath([code_loc '/src/library/']); 



% addpath('/Volumes/Data_CSNL/people/JWL/Toolbox/NIfTI_toolbox');
% addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/mi');
% addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Exp1_BOLD analysis/Behavior/'); % behavioral data
% addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/PopulationCoding_MT/Core analysis codes'); 
% addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Exp1_BOLD analysis/Library'); 
% addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Code Library'); 
% addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Tuned_fin_Nov'); 
% addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/PopulationCoding_MT'); 
% addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/BayesianSampling/library/'); % parameter setting
% addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7');

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







% Load behavioral data (model regressors)
fprintf('Load behavioral data.... \n'); 
ParamSetting_fmri;
% load([code_loc '/data/behavior/model_regressors.mat'], 'matSFM','matMimic','matSFM_PredCo');
load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7/MTanalysis_deep/LeeLee_PriorEffect_workspace_PredCo20171127.mat','matSFM','matMimic','matSFM_PredCo'); 
load('mDataMat_all'); 

ValSub =  1:21; 

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




locationROI = '/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/PopulationCoding_MT/ScanData_MTlocalization';
locationData = '/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/PopulationCoding_MT/ScanData_Main';
locationGLM = '/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/GLM_data_SFM';
locationSPM = '/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/GLM_finalVersion/GLM_PredictiveCodings/Model_PredCo_Weil/'; 
locationSPM_mimic = '/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/GLM_finalVersion/GLM_mimic/Model_U/'; 
locationRoot = '/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/PopulationCoding_MT';
locationSim = '/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/BayesianSampling/SimulationResult_3D/'; 

CalculateWhat = 1; 



for iSub = ValSub
    %% [2] Making BOLD time-series (population activity)
	fprintf(['Currently processing #' num2str(iSub) ', ' session.ID{iSub} '... \n']); 
    
    cd([locationData '/' session.ID{iSub} '/']);
    clear filename; 
    if iSub < 15
        filename.sfm{1} = ['s' session.ID{iSub} '_wr04SfM1.nii'];
        filename.sfm{2} = ['s' session.ID{iSub} '_wr05SfM2.nii'];
        filename.sfm{3} = ['s' session.ID{iSub} '_wr06SfM3.nii'];
        filename.sfm{4} = ['s' session.ID{iSub} '_wr07SfM4.nii'];
        
        filename.mimic{1} = ['s' session.ID{iSub} '_wr08SfMmimic.nii'];
    else
        filename.sfm{1} = ['s' session.ID{iSub} '_wr2SfM1.nii'];
        filename.sfm{2} = ['s' session.ID{iSub} '_wr4SfM2.nii'];
        filename.sfm{3} = ['s' session.ID{iSub} '_wr6SfM3.nii'];
        filename.sfm{4} = ['s' session.ID{iSub} '_wr8SfM4.nii'];

        filename.mimic{1} = ['s' session.ID{iSub} '_wr3Mimic1.nii'];
        filename.mimic{2} = ['s' session.ID{iSub} '_wr5Mimic2.nii'];
        filename.mimic{3} = ['s' session.ID{iSub} '_wr7Mimic3.nii'];
        filename.mimic{4} = ['s' session.ID{iSub} '_wr9Mimic4.nii'];
    end

    clear NiftiFile_sfm NiftiFile_mimic 
	% Loading Nifti files 
    for iS = 1:length(filename.sfm)
        NiftiFile_sfm{iS} = load_nii(filename.sfm{iS});
    end

    for iS = 1:length(filename.mimic)
        NiftiFile_mimic{iS} = load_nii(filename.mimic{iS});
    end
    [nx,ny,nz,nt] = size(NiftiFile_sfm{1}.img);
    Nifti_AlliS = zeros(nx, ny, nz, (nt)*4); 
    for iS = 1:4
        Nifti_AlliS(:,:,:,(1:nt) + (iS-1)*(nt)) = single(squeeze(NiftiFile_sfm{iS}.img(:,:,:,1:end)))-repmat(mean(NiftiFile_sfm{iS}.img,4),1,1,1,125);
    end
    
    Nifti_AlliS_rep = zeros(nx, ny, nz, (nt)*length(filename.mimic)); 
    for iS = 1:length(filename.mimic)
        Nifti_AlliS_rep(:,:,:,(1:nt) + (iS-1)*(nt)) = single(squeeze(NiftiFile_mimic{iS}.img(:,:,:,1:end)))-repmat(mean(NiftiFile_mimic{iS}.img,4),1,1,1,125);
    end
    
    % Preprocessing (HPF & percent signal)
    clear Nifti_AlliS_filt Nifti_AlliS_rep_filt; % reset
    Nifti_AlliS_filt = zeros(length(GrayCoord),500); 
	if iSub > 14
		Nifti_AlliS_rep_filt = zeros(length(GrayCoord),500); 
	else
		Nifti_AlliS_rep_filt = zeros(length(GrayCoord),125); 
	end
	
    for iVox = 1:length(GrayCoord)
        temp123 = squeeze(Nifti_AlliS(GrayCoord(iVox,1), GrayCoord(iVox,2), GrayCoord(iVox,3), :));
        temp124 = squeeze(Nifti_AlliS_rep(GrayCoord(iVox,1), GrayCoord(iVox,2), GrayCoord(iVox,3), :));
        if sum(abs(diff(temp123)))~=0
            % z-scoring
            Nifti_AlliS_filt(iVox,:) = ((temp123-nanmean(temp123))/nanstd(temp123))';
            Nifti_AlliS_rep_filt(iVox,:) = ((temp124-nanmean(temp124))/nanstd(temp124))';
        else
            Nifti_AlliS_filt(iVox,:) = nan(size(temp123))';
            Nifti_AlliS_rep_filt(iVox,:) = nan(size(temp124))';
        end
    end
    
    % Subtract run-average fluctuation
    n_1 = length(Nifti_AlliS_filt(1,:))/125;
    n_2 = length(Nifti_AlliS_rep_filt(1,:))/125;
    Nifti_AlliS_filt = Nifti_AlliS_filt - repmat(nanmean(reshape([Nifti_AlliS_filt Nifti_AlliS_rep_filt], length([Nifti_AlliS_filt Nifti_AlliS_rep_filt]), 125,(n_1+n_2)),3),1,n_1); 
    Nifti_AlliS_rep_filt = Nifti_AlliS_rep_filt - repmat(nanmean(reshape([Nifti_AlliS_filt Nifti_AlliS_rep_filt], length([Nifti_AlliS_filt Nifti_AlliS_rep_filt]), 125,(n_1+n_2)),3),1,n_2); 
    
    
    % Visualize the time-series of an exempler voxel. 
    set(figure(1),'position',[2 289 1040 516]); clf; 
    VoxRep = 1; VoxRep2 = 200; 
    SP = subplot('position',[0.08 0.85 0.2 0.1]); cla; hold on;
    plot(time_fmri, Nifti_AlliS_filt(VoxRep,1:125),'color','k','linewidth',1.5);
    plot(time_fmri, Nifti_AlliS_filt(VoxRep2,1:125),'color',[0 0 0]+0.6,'linewidth',1.5);
    xlabel('Time (TR)'); ylabel('BOLD (z)'); xlim([0 300]); ylim([-5 5]); 
    
    SP = subplot('position',[0.08 0.6 0.2 0.1]); cla; hold on;
    plot(time_fmri, PE_all{iSub}(1,1:125),'color','k'); xlim([0 300]); ylim([-0.1 1.1]); 
    ylabel('Expectation');
    
    SP = subplot('position',[0.08 0.45 0.2 0.1]); cla; hold on;
    plot(time_fmri, Percept_all{iSub}(1,1:125),'color','k'); ylim([-1.1 1.1]); 
    ylabel('Choice');
    
    SP = subplot('position',[0.08 0.3 0.2 0.1]); cla; hold on;
    plot(time_fmri, matSFM{1,iSub}.result_rs.UU,'color','k');
    ylabel('Surprise (U)'); ylim([0 1]); 
    
    SP = subplot('position',[0.08 0.15 0.2 0.1]); cla; hold on;
    plot(time_fmri, matSFM{1,iSub}.result_rs.EU,'color','k');
    xlabel('Time (TR)'); ylabel('Entropy (H)'); ylim([0 1]);  
    
    % save BOLD mat
    save([code_loc '/data/wholebrain/wb_bold_sub' num2str(iSub) '.mat'],'Nifti_AlliS_filt', 'Nifti_AlliS_rep_filt')

    SimpleCorr_sfm = nan(1,length(GrayCoord)); 
    SimplePval_sfm = nan(1,length(GrayCoord)); 
    PartCorr_sfm = nan(1,length(GrayCoord)); 
    PartPval_sfm = nan(1,length(GrayCoord)); 
    SimpleCorr_mimic = nan(1,length(GrayCoord)); 
    SimplePval_mimic = nan(1,length(GrayCoord)); 
    PartCorr_mimic = nan(1,length(GrayCoord)); 
    PartPval_mimic = nan(1,length(GrayCoord)); 
    for iVox = 1:length(GrayCoord)
        nanInd = ~isnan(Nifti_AlliS_filt(iVox,:) + PE_all{iSub} + Percept_all{iSub}); 
        if sum(nanInd) > 5
            [SimpleCorr_sfm(iVox), SimplePval_sfm(iVox)] = corr(Nifti_AlliS_filt(iVox,nanInd)', PE_all{iSub}(nanInd)'); 
            [PartCorr_sfm(iVox), temp] = partcorr(Nifti_AlliS_filt(iVox,nanInd)', PE_all{iSub}(nanInd)', Percept_all{iSub}(nanInd)');
            PartPval_sfm(iVox) = temp.p; 
        end
        nanInd = ~isnan(Nifti_AlliS_rep_filt(iVox,:) + PE_all_mimic{iSub} + Percept_all_mimic{iSub}); 
        if sum(nanInd) > 5
            [SimpleCorr_mimic(iVox), SimplePval_mimic(iVox)] = corr(Nifti_AlliS_rep_filt(iVox,nanInd)', PE_all_mimic{iSub}(nanInd)'); 
            [PartCorr_mimic(iVox),temp] = partcorr(Nifti_AlliS_rep_filt(iVox,nanInd)', PE_all_mimic{iSub}(nanInd)', Percept_all_mimic{iSub}(nanInd)');
            PartPval_mimic(iVox) = temp.p; 
        end
    end
    Corr_simple.sfm_corr = SimpleCorr_sfm;
    Corr_simple.sfm_pval = SimplePval_sfm;
    Corr_simple.mimic_corr = SimpleCorr_mimic;
    Corr_simple.mimic_pval = SimplePval_mimic;
    
    Corr_part.sfm_corr = PartCorr_sfm;
    Corr_part.sfm_pval = PartPval_sfm;
    Corr_part.mimic_corr = PartCorr_mimic;
    Corr_part.mimic_pval = PartPval_mimic;
    save([code_loc '/data/wholebrain/corr_sub' num2str(iSub) '.mat'],'Corr_simple', 'Corr_part')
end




for iSub = ValSub
    iSub
    load([code_loc '/data/wholebrain/wb_bold_sub' num2str(iSub) '.mat'],'Nifti_AlliS_filt', 'Nifti_AlliS_rep_filt')

    SimpleCorr_sfm = nan(1,length(GrayCoord)); 
    SimplePval_sfm = nan(1,length(GrayCoord)); 
    PartCorr_sfm = nan(1,length(GrayCoord)); 
    PartPval_sfm = nan(1,length(GrayCoord)); 
    SimpleCorr_mimic = nan(1,length(GrayCoord)); 
    SimplePval_mimic = nan(1,length(GrayCoord)); 
    PartCorr_mimic = nan(1,length(GrayCoord)); 
    PartPval_mimic = nan(1,length(GrayCoord)); 
    for iVox = 1:length(GrayCoord)
        nanInd = ~isnan(Nifti_AlliS_filt(iVox,:) + DV_all{iSub} + Percept_all{iSub}); 
        if sum(nanInd) > 5
            [SimpleCorr_sfm(iVox), SimplePval_sfm(iVox)] = corr(Nifti_AlliS_filt(iVox,nanInd)', DV_all{iSub}(nanInd)'); 
            [PartCorr_sfm(iVox), temp] = partcorr(Nifti_AlliS_filt(iVox,nanInd)', DV_all{iSub}(nanInd)', Percept_all{iSub}(nanInd)');
            PartPval_sfm(iVox) = temp.p; 
        end
%         nanInd = ~isnan(Nifti_AlliS_rep_filt(iVox,:) + PE_all_mimic{iSub} + Percept_all_mimic{iSub}); 
%         if sum(nanInd) > 5
%             [SimpleCorr_mimic(iVox), SimplePval_mimic(iVox)] = corr(Nifti_AlliS_rep_filt(iVox,nanInd)', PE_all_mimic{iSub}(nanInd)'); 
%             [PartCorr_mimic(iVox),temp] = partcorr(Nifti_AlliS_rep_filt(iVox,nanInd)', PE_all_mimic{iSub}(nanInd)', Percept_all_mimic{iSub}(nanInd)');
%             PartPval_mimic(iVox) = temp.p; 
%         end
    end
    save([code_loc '/data/wholebrain/DV/corr_sub' num2str(iSub) '.mat'],'Corr_simple', 'Corr_part')
    
    SimpleCorr_sfm = nan(1,length(GrayCoord)); 
    SimplePval_sfm = nan(1,length(GrayCoord)); 
    PartCorr_sfm = nan(1,length(GrayCoord)); 
    PartPval_sfm = nan(1,length(GrayCoord)); 
    SimpleCorr_mimic = nan(1,length(GrayCoord)); 
    SimplePval_mimic = nan(1,length(GrayCoord)); 
    PartCorr_mimic = nan(1,length(GrayCoord)); 
    PartPval_mimic = nan(1,length(GrayCoord)); 
    for iVox = 1:length(GrayCoord)
        nanInd = ~isnan(Nifti_AlliS_filt(iVox,:) + DV_all{iSub} + Percept_all{iSub}); 
        if sum(nanInd) > 5
            [SimpleCorr_sfm(iVox), SimplePval_sfm(iVox)] = corr(Nifti_AlliS_filt(iVox,nanInd)', Percept_all{iSub}(nanInd)'); 
%             [PartCorr_sfm(iVox), temp] = partcorr(Nifti_AlliS_filt(iVox,nanInd)', DV_all{iSub}(nanInd)', Percept_all{iSub}(nanInd)');
%             PartPval_sfm(iVox) = temp.p; 
        end
%         nanInd = ~isnan(Nifti_AlliS_rep_filt(iVox,:) + PE_all_mimic{iSub} + Percept_all_mimic{iSub}); 
%         if sum(nanInd) > 5
%             [SimpleCorr_mimic(iVox), SimplePval_mimic(iVox)] = corr(Nifti_AlliS_rep_filt(iVox,nanInd)', PE_all_mimic{iSub}(nanInd)'); 
%             [PartCorr_mimic(iVox),temp] = partcorr(Nifti_AlliS_rep_filt(iVox,nanInd)', PE_all_mimic{iSub}(nanInd)', Percept_all_mimic{iSub}(nanInd)');
%             PartPval_mimic(iVox) = temp.p; 
%         end
    end
    save([code_loc '/data/wholebrain/Decision/corr_sub' num2str(iSub) '.mat'],'Corr_simple', 'Corr_part')
end


%% Figure 1 
Fig_omfc_1









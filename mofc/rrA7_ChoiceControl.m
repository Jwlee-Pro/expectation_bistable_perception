
clear all; 

%% [1] Preparation 
    global gParam time_fmri session

% Subjects to analyze
    ValSub =  1:21; 

    
% General parameters 
    gParam.TR = 2.4;                % TR
    gParam.HD = 2;                  % Hemodynamic delay (for HD shift and population decoding) 
    gParam.EventWindow = 34;        % +-20sec window for Transition event 
    gParam.FigOn = 1;
	gParam.TrCrit = [2 250 2];     % Transition trial selection criterion (pre-stable, unstable, post-stable duration)
    time_fmri = (1:125)*gParam.TR;
    
    gParam.time_fmri_rs = time_fmri(1):0.1:time_fmri(end); 
    
    
    
    
% Data loading information & toolbox setting
fprintf('Data loading information & toolbox setting.... \n'); 
addpath('/Volumes/Data_CSNL/people/JWL/Toolbox/NIfTI_toolbox');
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/mi');
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Exp1_BOLD analysis/Behavior/'); % behavioral data
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/PopulationCoding_MT/Core analysis codes'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Exp1_BOLD analysis/Library'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Code Library'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Tuned_fin_Nov'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/PopulationCoding_MT'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/BayesianSampling/library/'); % parameter setting
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7');
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/GLMbased_PE'); 



ParamSetting;
cmap_x = copper(8); 
cmap_allmodel = [cmap_x(end-5:end,:); [255 111 207]/255 ; [154 206 81]/255]; 
cmap_choice=[102 204 255;204 204 204; 255 204 102]./255;
Common_fontsize = 10; 


% Load behavioral data 
fprintf('Load behavioral data.... \n'); 
load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7/MTanalysis_deep/LeeLee_PriorEffect_workspace_PredCo20171127.mat','matSFM','matMimic','matSFM_PredCo'); 

clear DataMat
for iSub = ValSub
    for iS = 1:4
        DataMat.Regressor.percept{iSub,iS}  = matSFM{iS,iSub}.result_rs.Percept'; 
        DataMat.Regressor.PE{iSub,iS}       = matSFM{iS,iSub}.result_rs.matPrior'; 
    end
end


for iSub = ValSub
    load(['/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/GLMbased_PE/FittedResult_WB_' num2str(iSub) '.mat'],'RegTS'); 
    
    DataMat.PEcorr_ChoiceControl = nan(length(RegTS{iSub,1}),4);
    DataMat.PEpval_ChoiceControl = nan(length(RegTS{iSub,1}),4);
    % Fit the time-series to the model  
    for iS = 1:4 
        fprintf(['Session #' num2str(iSub) ', iS = ' num2str(iS) '...\n'])
        for iVox = 1:length(RegTS{iSub,iS}(:,1)) 
            if sum(isnan(RegTS{iSub,iS}(iVox,:))) < 5
                targetTS = RegTS{iSub,iS}(iVox,:); 
                tempX = [targetTS; DataMat.Regressor.PE{iSub,iS}'; DataMat.Regressor.percept{iSub,iS}'];
                ValInd = ~isnan(mean(tempX)); 
                
                [DataMat.PEcorr_ChoiceControl(iVox,iS),a2] = partcorr(targetTS', DataMat.Regressor.PE{iSub,iS}, DataMat.Regressor.percept{iSub,iS}); 
                DataMat.PEpval_ChoiceControl(iVox,iS) = a2.p; 
            end
        end
    end
    save(['/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/GLMbased_PE/PartCorrelation_WB_' num2str(iSub) '.mat'],'DataMat')
end


% % Merge runs 
% for iSub = ValSub
%     load(['/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/GLMbased_PE/FittedResult_WB_' num2str(iSub) '.mat'],'RegTS'); 
%     
%     DataMat_runMerge.PEcorr_ChoiceControl = nan(length(RegTS{iSub,1}),1);
%     DataMat_runMerge.PEpval_ChoiceControl = nan(length(RegTS{iSub,1}),1);
%     
%     % Fit the time-series to the model 
%     fprintf(['Session #' num2str(iSub) '...\n'])
%     
%     models = []; 
%     for iS = 1:4
%         models = [models [DataMat.Regressor.PE{iSub,iS}'; DataMat.Regressor.percept{iSub,iS}']];
%     end
%     
%     for iVox = 1:length(RegTS{iSub,iS}(:,1)) 
%         targetTS = [RegTS{iSub,1}(iVox,:) RegTS{iSub,2}(iVox,:) RegTS{iSub,3}(iVox,:) RegTS{iSub,4}(iVox,:)]; 
%         ValInd = ~isnan(mean([models; targetTS])); 
%         if sum(ValInd) > 200
%             [DataMat_runMerge.PEcorr_ChoiceControl(iVox,1), a2] = partcorr(targetTS(ValInd)', models(1,ValInd)', models(2,ValInd)'); 
%             DataMat_runMerge.PEpval_ChoiceControl(iVox,1) = a2.p; 
%         end
%     end
%     save(['/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/GLMbased_PE/PartCorrelation_WB_' num2str(iSub) '_concat.mat'],'DataMat_runMerge')
% end

%% Permutation test 

for iSub = ValSub
    load(['/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/GLMbased_PE/FittedResult_WB_' num2str(iSub) '.mat'],'RegTS'); 
    axxx = [];
    for iS = 1:4
        axxx = [axxx nanmean(MainDatax{2,iS})]; 
    end
    plot(axxx);
    DataMat.PEcorr_ChoiceControl = nan(length(RegTS{iSub,1}),4);
    DataMat.PEpval_ChoiceControl = nan(length(RegTS{iSub,1}),4);
    % Fit the time-series to the model  
    for iS = 1:4 
        fprintf(['Session #' num2str(iSub) ', iS = ' num2str(iS) '...\n'])
        for iVox = 1:length(RegTS{iSub,iS}(:,1)) 
            if sum(isnan(RegTS{iSub,iS}(iVox,:))) < 5
                targetTS = RegTS{iSub,iS}(iVox,:); 
                tempX = [targetTS; DataMat.Regressor.PE{iSub,iS}'; DataMat.Regressor.percept{iSub,iS}'];
                ValInd = ~isnan(mean(tempX)); 
                
                [DataMat.PEcorr_ChoiceControl(iVox,iS),a2] = partcorr(targetTS', DataMat.Regressor.PE{iSub,iS}, DataMat.Regressor.percept{iSub,iS}); 
                DataMat.PEpval_ChoiceControl(iVox,iS) = a2.p; 
            end
        end
    end
    save(['/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/GLMbased_PE/PartCorrelation_WB_' num2str(iSub) '.mat'],'DataMat')
end

asdf

%% After calculating the correlation values 

for iSub = 1:21
%     load(['/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/GLMbased_PE/PartCorrelation_WB_' num2str(iSub) '_concat.mat'],'DataMat_runMerge'); 
%     DataMat_merge.PEpval_afReg{iSub} = DataMat_runMerge.PEpval_ChoiceControl; 
%     DataMat_merge.PEcorr_afReg{iSub} = DataMat_runMerge.PEcorr_ChoiceControl; 

    load(['PartCorrelation_WB_' num2str(iSub)],'DataMat');
    DataMat_merge.PEpval_afReg{iSub} = DataMat.PEpval_ChoiceControl; 
    DataMat_merge.PEcorr_afReg{iSub} = DataMat.PEcorr_ChoiceControl; 
end

% for iSub = 1:21
%     numvv(iSub) = sum(isnan(DataMat_merge.PEcorr_afReg{iSub})); 
% end


MatSig      = nan(length(DataMat_merge.PEpval_afReg), length(DataMat_merge.PEpval_afReg{1})); 
MatAvg_corr = nan(length(DataMat_merge.PEpval_afReg), length(DataMat_merge.PEpval_afReg{1})); 
MatAvg_var  = nan(length(DataMat_merge.PEpval_afReg), length(DataMat_merge.PEpval_afReg{1})); 
P_crit = 0.05; 

for iSub = 1:21
    MatSig(iSub,:) = (nanmean(DataMat_merge.PEpval_afReg{iSub} < P_crit,2) > 0.7); % At least 75%...
    MatAvg_corr(iSub,:) = nanmean(DataMat_merge.PEcorr_afReg{iSub},2); 
    MatAvg_var(iSub,:) = nanmean(DataMat_merge.PEcorr_afReg{iSub}.^2,2);     
end

% for iSub = 1:21
%     kV(iSub) = min(abs(DataMat_merge.PEcorr_afReg{iSub}((DataMat_merge.PEpval_afReg{iSub}(:,1) < P_crit),1))); 
% end


cd('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis'); 
ROI_Gray = load_nii('grey_dimMatched.nii');
[nx,ny,nz,nt] = size(ROI_Gray.img);
load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/GLMbased_PE/GrayCoordinates'); 

[h,p,ci,stats] = ttest(MatAvg_var); 
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p, 0.0001);
sum(h)
% min(adj_p(h==1))
% max(adj_p(h==1))


tt = mean(MatSig); 
tt1 = mean(MatAvg_corr); 
tt2 = mean(MatAvg_var); 

tempI = nan(nx,ny,nz); 
tempI1 = nan(nx,ny,nz); 
tempI2 = nan(nx,ny,nz); 
tempIx = nan(nx,ny,nz); 
tempIx_crit = nan(nx,ny,nz); 
tempIx_crit_ff = nan(nx,ny,nz); 
for iVox = 1:120948
    tempI(GrayCoord(iVox,1), GrayCoord(iVox,2), GrayCoord(iVox,3)) = tt(iVox); 
    tempI1(GrayCoord(iVox,1), GrayCoord(iVox,2), GrayCoord(iVox,3)) = tt1(iVox); 
    tempI2(GrayCoord(iVox,1), GrayCoord(iVox,2), GrayCoord(iVox,3)) = tt2(iVox); 
    if tt2(iVox) > 0.02
        tempIx(GrayCoord(iVox,1), GrayCoord(iVox,2), GrayCoord(iVox,3)) = tt2(iVox); 
        if adj_p(iVox) < 0.0007
            tempIx_crit_ff(GrayCoord(iVox,1), GrayCoord(iVox,2), GrayCoord(iVox,3)) = tt2(iVox); 
        end
    end
    if adj_p(iVox) < 0.0007
        tempIx_crit(GrayCoord(iVox,1), GrayCoord(iVox,2), GrayCoord(iVox,3)) = tt2(iVox); 
    end
end

100*sum(tt2>(kV^2))/length(tt2)


% ROI_Gray.img = tempI;
% save_nii(ROI_Gray,'sessconcatChoiceC_SigVox_density.nii'); 
% 
% ROI_Gray.img = tempI1;
% save_nii(ROI_Gray,'sessconcatChoiceC_SigVox_corr.nii'); 
% 
% ROI_Gray.img = tempI2;
% save_nii(ROI_Gray,'sessconcatChoiceC_SigVox_var.nii'); 
% 
% ROI_Gray.img = tempIx;
% save_nii(ROI_Gray,'sessconcatChoiceC_SigVox_var_thr.nii'); 

% ROI_Gray.img = tempI;
% save_nii(ROI_Gray,'ChoiceC_SigVox_density.nii'); 
% 
% ROI_Gray.img = tempI1;
% save_nii(ROI_Gray,'ChoiceC_SigVox_corr.nii'); 
% 
% ROI_Gray.img = tempI2;
% save_nii(ROI_Gray,'ChoiceC_SigVox_var.nii'); 
% 
% ROI_Gray.img = tempIx;
% save_nii(ROI_Gray,'ChoiceC_SigVox_var_thr.nii'); 



ROI_Gray.img = tempI; 
save_nii(ROI_Gray,'ChoiceC_SigVox_density2023Jan.nii'); 

ROI_Gray.img = tempI1; 
save_nii(ROI_Gray,'ChoiceC_SigVox_corr2023Jan.nii'); 

ROI_Gray.img = tempI2; 
save_nii(ROI_Gray,'ChoiceC_SigVox_var2023Jan.nii'); 

ROI_Gray.img = tempIx; 
save_nii(ROI_Gray,'ChoiceC_SigVox_var_thr2023Jan.nii'); 

ROI_Gray.img = tempIx_crit; 
save_nii(ROI_Gray,'ChoiceC_SigVox_var_thr2023Jan_final.nii'); 


ROI_Gray.img = tempIx_crit_ff; 
save_nii(ROI_Gray,'ChoiceC_SigVox_var_thr2023Feb_final.nii'); 

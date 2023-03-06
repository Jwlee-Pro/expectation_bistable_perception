
% Data loading information & toolbox setting
fprintf('Data loading information & toolbox setting.... \n'); 
addpath('/Volumes/Data_CSNL/people/JWL/Toolbox/NIfTI_toolbox');
addpath(genpath('/Volumes/Data_CSNL/people/JWL/Toolbox/libsvm-3.21'));
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Exp1_BOLD analysis/Behavior/'); % behavioral data
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/PopulationCoding_MT/Core analysis codes'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Exp1_BOLD analysis/Library'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Code Library'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Tuned_fin_Nov'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/PopulationCoding_MT'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/BayesianSampling/library/'); % parameter setting
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7');
 
 
ParamSetting; 
cmap_x = copper(8); 
cmap_allmodel = [cmap_x(end-5:end,:); [255 111 207]/255 ; [154 206 81]/255]; 
cmap_choice=[102 204 255;204 204 204; 255 204 102]./255;
cmap_choiceStrong = [226 180 90;150 150 150; 88 178 221]./255; 
Common_fontsize = 10; 



% General parameters 
    fprintf('General parameters setting.... \n'); 
    gParam.TR = 2.4;                % TR
    gParam.HD = 2;                  % Hemodynamic delay (for HD shift and population decoding) 
    gParam.EventWindow = 34;        % +-20sec window for Transition event 
    gParam.FigOn = 1;
	gParam.TrCrit = [2 250 2];     % Transition trial selection criterion (pre-stable, unstable, post-stable duration)
    time_fmri = (1:125)*gParam.TR;
    gParam.time_fmri_rs = time_fmri(1):0.1:time_fmri(end); 


% Load behavioral data 
fprintf('Load behavioral data.... \n'); 
ValSub =  1:21; 
load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7/MTanalysis_deep/LeeLee_PriorEffect_workspace_PredCo20171127.mat','matSFM','matMimic','matSFM_PredCo'); 

Percept_all = cell(1,21); Percept_all_mimic = cell(1,21); 
PE_all = cell(1,21); PE_all_mimic = cell(1,21); 
nSess = [ones(1,14) ones(1,7)*4];

for iSub = ValSub
    for iS = 1:4
        Percept_all{iSub} = [Percept_all{iSub} matSFM{iS,iSub}.result_rs.Percept]; 
        PE_all{iSub} = [PE_all{iSub} matSFM{iS,iSub}.result_rs.matPrior];         
    end
    
    for iS = 1:nSess(iSub)
        Percept_all_mimic{iSub} = [Percept_all_mimic{iSub} matMimic{iS,iSub}.result_rs.Percept]; 
        PE_all_mimic{iSub} = [PE_all_mimic{iSub} matMimic{iS,iSub}.result_rs.matPrior];         
    end
end


for iSub = 1:21
    nanInd = isnan(PE_all{iSub});
    if sum(nanInd) ~= 0
        PE_all{iSub}(nanInd) = PE_all{iSub}([nanInd(2:length(PE_all{iSub})) 0]==1);
    end
end
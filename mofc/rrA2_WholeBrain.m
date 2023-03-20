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


code_loc = '/Users/joonwon/Dropbox/Workspace_2022/Structure_from_motion/expectation_bistable_perception'; 


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
load([code_loc '/data/behavior/model_regressors.mat'], 'matSFM','matMimic','matSFM_PredCo');

ValSub =  1:21; 

Percept_all = cell(1,21); Percept_all_mimic = cell(1,21); 
PE_all = cell(1,21); PE_all_mimic = cell(1,21); 
U_all = cell(1,21); H_all = cell(1,21);
nSess = [ones(1,14) ones(1,7)*4];
for iSub = ValSub
    for iS = 1:4
        Percept_all{iSub} = [Percept_all{iSub} matSFM{iS,iSub}.result_rs.Percept]; 
        PE_all{iSub} = [PE_all{iSub} matSFM{iS,iSub}.result_rs.matPrior];         
        H_all{iSub} = [H_all{iSub} matSFM{iS,iSub}.result_rs.EU];         
        U_all{iSub} = [U_all{iSub} matSFM{iS,iSub}.result_rs.UU];         
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
    [d, hdr] = cbiReadNifti(filename.sfm{iS},{[],[],[],[]},'native');
    d = double(d);
    % High-Pass time series (used local mean)
    highpassPeriod = nTotalFrame / nTotalCycle;
    htSeries = removeLowFqDrift(tSeries, highpassPeriod);

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
    
    drawnow; 
    asdf
    
    clear MatFit; 
    if CalculateWhat == 1
        % Fit the time-series to the model 
        for iS = 1:4
            fprintf(['Session #' num2str(iSub) ', iS = ' num2str(iS) '...\n'])
            clear seedTS;
            seedTS(1,:) = matSFM{iS,iSub}.result_rs.matPrior;   % PE
            seedTS(2,:) = matSFM{iS,iSub}.result_rs.Percept;    % Choice
            seedTS(3,:) = matSFM{iS,iSub}.result_rs.UU;         % U
            seedTS(4,:) = 1:125;                                % Trend
            
            for iVox = 1:length(Nifti_AlliS_filt(:,1)) 
                if sum(isnan(Nifti_AlliS_filt(iVox,:))) < 5
                    targetTS = Nifti_AlliS_filt(iVox,(1:125) + (125*(iS-1))); 
                    tempX = [seedTS; targetTS];
                    ValInd = ~isnan(mean(tempX)); 

                    paramModel1 = [0 0 1 0 0]; 
                    [paramout_temp, fval_temp] = fminsearchbnd( @(param) fitModel_PE(param, targetTS(ValInd), seedTS(:,ValInd)), ...
                        paramModel1, [-10 -10 0 -10 -10] ,[10 10 10 10 10]);
                    MatFit{iS}(iVox,:) = paramout_temp; 
                else
                    MatFit{iS}(iVox,:) = nan(1,5); 
                end
            end
        end
    save(['/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/FittedResult_wholeBrain_' num2str(iSub) '.mat'],'MatFit')

%     % Calculate the correlation
%     MatCorr_sfm = nan(1,length(GrayCoord)); 
%     MatCorr_mimic = nan(1,length(GrayCoord)); 
%     
    elseif CalculateWhat == 2

        for iVox = 1:length(GrayCoord)
            nanInd = ~isnan(Nifti_AlliS_filt(iVox,:) + PE_all{iSub} + Percept_all{iSub}); 
            if sum(nanInd) > 5
                [SimpleCorr_sfm(iVox), SimplePval_sfm(iVox)] = corr(Nifti_AlliS_filt(iVox,nanInd)', PE_all{iSub}(nanInd)'); 
                MI(iVox) = calculate_mi(Nifti_AlliS_filt(iVox,nanInd),PE_all{iSub}(nanInd))
%                 [MatCorr_sfm(iVox), temp] = partcorr(Nifti_AlliS_filt(iVox,nanInd)', PE_all{iSub}(nanInd)', Percept_all{iSub}(nanInd)');
%                 MatPval_sfm(iVox) = temp.p; 
            end
            nanInd = ~isnan(Nifti_AlliS_rep_filt(iVox,:) + PE_all_mimic{iSub} + Percept_all_mimic{iSub}); 
            if sum(nanInd) > 5
                [SimpleCorr_mimic(iVox), SimplePval_mimic(iVox)] = corr(Nifti_AlliS_rep_filt(iVox,nanInd)', PE_all_mimic{iSub}(nanInd)'); 
%                 [MatCorr_mimic(iVox),temp] = partcorr(Nifti_AlliS_rep_filt(iVox,nanInd)', PE_all_mimic{iSub}(nanInd)', Percept_all_mimic{iSub}(nanInd)');
%                 MatPval_mimic(iVox) = temp.p; 
            end
        end
        MatCorr{iSub}.sfm_corr = SimpleCorr_sfm;
        MatCorr{iSub}.sfm_pval = SimplePval_sfm;
        
        MatCorr{iSub}.mimic_corr = SimpleCorr_mimic;
        MatCorr{iSub}.mimic_pval = SimplePval_mimic;
        
    end
%     
%     
%     % Save Data
%     save(['/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/CorrelationMaps/AllCorrs' num2str(iSub) '_smooth.mat'],'MatCorr_sfm', 'MatCorr_mimic','SimpleCorr_sfm', 'SimpleCorr_mimic','MatPval_sfm','MatPval_mimic','SimplePval_sfm','SimplePval_mimic');
% %     save(['/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/FinalAnalysis/TS_OMPFC' num2str(iSub) '.mat'],'Nifti_AlliS_filt','Nifti_AlliS_rep_filt');
end




%% Analyze the fitting result 
SubRep = 4; 






% 
% 
% SP = subplot('position',[0.1300 0.75 0.1237 0.2]); cla; hold on;
% % Explained variance 
% SubRep = 4; 
% p_crit = 0.01; 
% xaxis_a = -1:0.02:1; 
% a = hist(MatCorr{SubRep}.sfm_corr,xaxis_a);
% plot(xaxis_a, a,'color','k'); 
% b = hist(MatCorr{SubRep}.sfm_corr(MatCorr{SubRep}.sfm_pval < p_crit),xaxis_a);
% plot(xaxis_a, b,'color','k','linewidth',1); 
% V = [xaxis_a; b]'; F = 1:length(V); 
% patch('Faces',F,'Vertices',V); 
% xlim([-1 1]); 
% xlabel('corr(PE, BOLD)'); ylabel('Number of voxels'); 
% set(SP, 'box', 'off', 'TickDir','out','XTick', [-1 0 1],'FontSize', 10,'ycolor','k','xcolor','k')
% 
% 
% MaxInd = find(MatCorr{SubRep}.sfm_corr == max(MatCorr{SubRep}.sfm_corr)); 
% MinInd = find(MatCorr{SubRep}.sfm_corr == min(MatCorr{SubRep}.sfm_corr)); 
% figure(2); clf; 
% subplot(2,4,1:3); plot(PE_all{SubRep},'k');
% 
% subplot(2,4,5:7); plot(Nifti_AlliS_filt(MaxInd, :),'r'); hold on; 
%  plot(Nifti_AlliS_filt(MinInd, :),'b'); 
% subplot(3,4,[4 8 12]); plot(PE_all{SubRep}, Nifti_AlliS_filt(MaxInd, :),'r.'); hold on; 
% plot(PE_all{SubRep}, Nifti_AlliS_filt(MinInd, :),'b.'); 


















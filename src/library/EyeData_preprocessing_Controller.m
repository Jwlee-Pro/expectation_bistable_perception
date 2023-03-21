%% Preprocessing 
clear all; 

addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Code Library/'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/New_Version_analysis_AllSbj_20151022/5.Pupil Analysis/Eyedata'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Exp2_Pupil analysis/Behavior'); 

% Scan data for main and mimic runs (1 run each) : Behavior 
scan.sessionID{1} = 'XJWL001';
scan.sessionID{2} = 'XJWL002';
scan.sessionID{3} = 'XKWC045';
scan.sessionID{4} = 'XJWL005';
scan.sessionID{5} = 'XJWL006';
scan.sessionID{6} = 'XJWL007';
scan.sessionID{7} = 'XJWL008';
scan.sessionID{8} = 'XJWL009';
scan.sessionID{9} = 'XKWC068';

% savePath = '/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Exp2_Pupil analysis/Eye_data_processed_20160626/'; 
savePath = '/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Exp2_Pupil analysis/Eye_data_processed_20191029/'; 

for iSess = 1:length(scan.sessionID)
    for iRun = 1:12
        if (exist([scan.sessionID{iSess} '_scanData_' num2str(iRun) '.mat'], 'file') == 2) && (exist(['eyetrack_' scan.sessionID{iSess} '_' num2str(iRun) 's.mat'], 'file') == 2) && (exist(['eyetrack_' scan.sessionID{iSess} '_' num2str(iRun) 'm.mat'], 'file') == 2)
            load([scan.sessionID{iSess} '_scanData_' num2str(iRun) '.mat'], 'paramSFM','paramMimic'); 
            % Eye data for main and mimic : Eye
            Eyedata{1}= load(['eyetrack_' scan.sessionID{iSess} '_' num2str(iRun) 's.mat']); % main run
            Eyedata{2}= load(['eyetrack_' scan.sessionID{iSess} '_' num2str(iRun) 'm.mat']); % mimic run

            corrected   = EyeData_preprocessing_Fin(Eyedata{1}.Eyelink_data, paramSFM, iSess);
            corrected_m = EyeData_preprocessing_Fin(Eyedata{2}.Eyelink_data, paramMimic, iSess);
        
            save([savePath, 'eyeData_processed_' scan.sessionID{iSess} '_' num2str(iRun) 's'],'corrected'); 
            save([savePath, 'eyeData_processed_' scan.sessionID{iSess} '_' num2str(iRun) 'm'],'corrected_m'); 

        end
    end
end





            
            
            
            
            
            
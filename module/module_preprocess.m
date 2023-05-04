%% Preprocessing module (overview)
% 1. Standardize the BOLD signal 
% 2. Slow fMRI fluctuation correction

step_HPF = 1; 
step_CSF = 1; 
step_WM  = 1; 
step_meanRegress = 1; 



%% Standardize the BOLD signal 
% Run-wise mean centering 
% Normalize with std 


%% Slow fMRI fluctuation correction
% - This is not related to the neural response
% - Potential cause: MRI temporature, fatigue, ...

% Non-neural activity: outside of brain / csf 
% Regress out, subtraction 


%% Main code 
% fMRI data is head motion-corrected, time-corrected, warped/resliced to MNI152 standard space
locationData = '/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/PopulationCoding_MT/ScanData_Main';

for iSub = 1:length(session.ID)
	fprintf(['Preprocessing #' num2str(iSub) ', ' session.ID{iSub} ' out of' num2str(length(session.ID)) ' ... \n']); 
    
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

    % mean correction 
    Nifti_AlliS = zeros(nx, ny, nz, (nt)*4); 
    for iS = 1:4
        Nifti_AlliS(:,:,:,(1:nt) + (iS-1)*(nt)) = single(squeeze(NiftiFile_sfm{iS}.img(:,:,:,1:end)))-repmat(mean(NiftiFile_sfm{iS}.img,4),1,1,1,125);
    end
    
    Nifti_AlliS_rep = zeros(nx, ny, nz, (nt)*length(filename.mimic)); 
    for iS = 1:length(filename.mimic)
        Nifti_AlliS_rep(:,:,:,(1:nt) + (iS-1)*(nt)) = single(squeeze(NiftiFile_mimic{iS}.img(:,:,:,1:end)))-repmat(mean(NiftiFile_mimic{iS}.img,4),1,1,1,125);
    end
    
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
    
    % save BOLD mat
    save([code_loc '/data/wholebrain/wb_bold_sub' num2str(iSub) '.mat'],'Nifti_AlliS_filt', 'Nifti_AlliS_rep_filt')

    
end

































function [] = glme1_NuisanceOut(serverpath,Conds)
%
subjID{01} = 'ELJ110405';
subjID{02} = 'MK101110';
subjID{03} = 'XS13_101110';
subjID{04} = 'XS14_101101';
subjID{05} = 'XS15_101103';
subjID{06} = 'XS16_101105';
subjID{07} = 'XS17_101112';
subjID{08} = 'XS18_101103';
subjID{09} = 'XS20_101105';
subjID{10} = 'XS21_101102';
subjID{11} = 'XS42_110325';
subjID{12} = 'XS43_110324';
subjID{13} = 'XS44_110329';
subjID{14} = 'XS45_110328';
subjID{15} = 'XS46_110323';
subjID{16} = 'XS47_110401';
subjID{17} = 'XS49_110331';
subjID{18} = 'YMS101122';

nROI = Conds.nROI;
ROIname = cell(nROI,1);
for iROI = 1:nROI
    ROIname{iROI} = num2str(iROI);
end
%
if Conds.CSFout == 1
    nuisance = cell(18,1);
    for iSub = 1:18
        load(['/Volumes/Data_CS' serverpath '/people/LHS/PDM whole brain/Dartel/CSFsignal/' subjID{iSub} '.mat'])
        nuisance{iSub} = nuibold';
    end
    nuisancesave = '';
else
    nuisancesave = '_no';
end
%
cutoff = Conds.cutoff;
icutoff = cutoff*ones(18,8);
hpfsave = ['hpf' num2str(cutoff) 's_'];
%
if isempty(dir(['/Volumes/Data_CS' serverpath '/people/LHS/results_ver11/BOLD/glme1_NuisanceOut/'])) == 1
    mkdir(['/Volumes/Data_CS' serverpath '/people/LHS/results_ver11/BOLD/glme1_NuisanceOut/'])
end
%%
cd(['/Volumes/Data_CS' serverpath '/people/LHS/results_ver11/BOLD/glme1_NuisanceOut/'])
for iROI = 1:nROI
    isexist1 = exist(['/Volumes/Data_CS' serverpath '/people/LHS/results_ver11/BOLD/glme1_NuisanceOut/' hpfsave ROIname{iROI} nuisancesave '_computing.mat'],'file');
    isexist2 = exist(['/Volumes/Data_CS' serverpath '/people/LHS/results_ver11/BOLD/glme1_NuisanceOut/' hpfsave ROIname{iROI} nuisancesave '.mat'],'file');
    %     							isexist1 = 0; isexist2 = 0;
    if isexist1 == 0 && isexist2 ==0
        save(['./' hpfsave ROIname{iROI} nuisancesave '_computing.mat'],'iROI')
        allBOLD = cell(18,1);
        for iSub = 1:18
            disp(['iROI = ' ROIname{iROI} ', BOLD loading, iSub = ' num2str(iSub) '...'])
            % load BOLD data
            load(['/Volumes/Data_CS' serverpath '/people/LHS/PDM whole brain/Dartel/functional/Parcellation/nROI' num2str(nROI) '/S' num2str(iSub) '_' num2str(iROI) '.mat'])
            allBOLD{iSub} = maskedBOLD2d{1};
        end
        % detrend -> high pass filter -> percent signal -> CFS out
        for iSub = 1:18
            for iRun = 1:8
                % 1. detrend
                iT		= (6*(iRun-1)*26+7):6*iRun*26;
                ibold	= allBOLD{iSub}(:,iT)';
                mbold	= repmat(mean(ibold),150,1);
                ibold	= detrend(ibold);
                if Conds.CSFout == 1
                    inui	= nuisance{iSub}(iT,:);
                    mnui	= repmat(mean(inui),150,1);
                    inui	= detrend(inui);
                end
                % 2. high-pass filter
                if sum(isnan(cutoff)) == 0
                    Fs=1/2.2; n=5; Fn=Fs/2; ftype='high';
                    Wn=1/icutoff(iSub,iRun);
                    [b,a]	= butter(n,Wn/Fn,ftype);
                    ibold	 = filtfilt(b,a,ibold);
                    if Conds.CSFout == 1
                        inui	= filtfilt(b,a,inui);
                    end
                end
                % 3. percent signal
                ibold = ibold./mbold*100;
                iVoxs = find(isnan(sum(ibold))==0);
                if Conds.CSFout == 1
                    % 4. nuisance out
                    inui = inui./mnui*100;
                    inui = mean(inui,2);
                    for iVox = iVoxs
                        jbold	= ibold(:,iVox);
                        if isnan(sum(jbold)) == 0
                            jbold	= spm_orth([inui jbold]);
                            jbold	= jbold(:,end);
                            ibold(:,iVox) = jbold;
                        end
                    end
                end
                allBOLD{iSub}(:,iT) = ibold';
            end
            fprintf('phase check... cutoff=%d, iROI=%1d, iSub=%d\n',cutoff,iROI,iSub)
        end
        maskedBOLD2d{1} = allBOLD;
        save(['./' hpfsave ROIname{iROI} nuisancesave '.mat'],'maskedBOLD2d')
        delete(['./' hpfsave ROIname{iROI} nuisancesave '_computing.mat'])
    end
end
end

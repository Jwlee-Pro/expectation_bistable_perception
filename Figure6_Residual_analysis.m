%% Residual analysis 


clear all; 
addpath('/Users/joonwon/Desktop/Matlab_local/Brain_residualAnalysis_2018June/PleaseWork/NeatVersion_library'); 
addpath('/Users/joonwon/Desktop/Matlab_local/Brain_residualAnalysis_2018June/PleaseWork'); 



%% Load behavioral data 
Setting_param_bhv; 
load('TS_MT_all.mat');
load('TS_OMPFC_all.mat');
load('TS_FP_all15.mat');

Generating_model_TS; 





%% Quantify the residual expectation (PE at OMPFC)

    PE_reg = cell(1,21); Dec_reg = cell(1,21); 
    for iSub = 1:21
        for iS = 1:4
            PE_reg{iSub} = [PE_reg{iSub} mDataMat.PE_conv{iSub,iS}]; 
            Dec_reg{iSub} = [Dec_reg{iSub} mDataMat.Decision_conv{iSub,iS}];         
        end
        CorrVal.OMPFC_PE{iSub} = corr(DataMat_OMPFC{iSub}.sfm', PE_reg{iSub}');
        CorrVal.MT_D{iSub} = corr(DataMat_MT{iSub}.sfm', Dec_reg{iSub}');
    end

    % Regressing out the U and H 
    for iSub = 1:21
        for iS = 1:4
            for iVox = 1:length(DataMat_OMPFC{iSub}.sfm(:,1))
                TStemp = DataMat_OMPFC{iSub}.sfm(iVox,(1:125) + (iS-1)*125); 
                [b,bint,r,rint,stats] = regress(TStemp', [mDataMat.U_conv{iSub,iS}' ones(size(mDataMat.U_conv{iSub,iS}'))]); 
                DataMat_OMPFC{iSub}.sfm_reg(iVox,(1:125) + (iS-1)*125) = r';
            end
        end
    end

    % Normalizing each bin's time-series
    for iSub = 1:21
        DataMat_OMPFC{iSub}.sfm_norm = (DataMat_OMPFC{iSub}.sfm_reg - repmat(mean(DataMat_OMPFC{iSub}.sfm_reg,2), 1,500))./repmat(std(DataMat_OMPFC{iSub}.sfm_reg,0,2), 1,500);
        DataMat_MT{iSub}.sfm_norm = (DataMat_MT{iSub}.sfm - repmat(mean(DataMat_MT{iSub}.sfm,2), 1,500))./repmat(std(DataMat_MT{iSub}.sfm,0,2), 1,500);
    end

    % Dissociate choice-explained variance and PE-explained variance
    clear partCorrMat partCorrMat_mt
    for iSub = 1:21
        for iVox = 1:length(DataMat_OMPFC{iSub}.sfm(:,1))
            partCorrMat.D(iSub,iVox) = partcorr(DataMat_OMPFC{iSub}.sfm_reg(iVox,:)', Dec_reg{iSub}', PE_reg{iSub}'); 
            partCorrMat.PE(iSub,iVox) = partcorr(DataMat_OMPFC{iSub}.sfm_reg(iVox,:)', PE_reg{iSub}', Dec_reg{iSub}'); 
        end

        signPE = (CorrVal.OMPFC_PE{iSub}>0)*2 -1; 
        partCorrMat.PE_signed(iSub,:) = partCorrMat.PE(iSub,:).*signPE'; 
        partCorrMat.D_signed(iSub,:) = partCorrMat.D(iSub,:).*signPE'; 
    end
    [aa1,~]=ttest((partCorrMat.D_signed));
    [aa2,~]=ttest((partCorrMat.PE_signed));
    for iVox = 1:length(aa1)
        if (aa1(iVox) == 1) && (aa2(iVox) == 1)
            indK(iVox) = 4; 
        elseif (aa1(iVox) == 1) && (aa2(iVox) ~= 1)
            indK(iVox) = 2; 
        elseif (aa1(iVox) ~= 1) && (aa2(iVox) == 1)
            indK(iVox) = 3; 
        else
            indK(iVox) = 1; % Both insignificant
        end
    end
    
    % Calculating the slope 
    for iSub = 1:21
        for iT = 1:length(DataMat_OMPFC{iSub}.sfm_norm(1,:))
            tp = polyfit(CorrVal.OMPFC_PE{iSub}(indK==3), DataMat_OMPFC{iSub}.sfm_norm((indK==3),iT), 1); 
            slopeDec{iSub}(1,iT) = tp(1); 
            tp = polyfit(CorrVal.OMPFC_PE{iSub}, DataMat_OMPFC{iSub}.sfm_norm(:,iT), 1); 
            stemp{iSub}(1,iT) = tp(1); 
        end
    end
    
    % MT decoded data
    load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7/MTanalysis_deep/LeeLee_PriorEffect_workspace_PredCo_2018June.mat'); 

    % Calculate the residual expectation
    MTorg = cell(1,21); 
    for iSub = 1:21
        for iS = 1:4
            MTorg{iSub} = [MTorg{iSub} dec_final_fin{4}{iSub,iS}]; 
        end
        [b,bint,resiDec{iSub},rint,stats] = regress(slopeDec{iSub}', [PE_reg{iSub}' ones(size(PE_reg{iSub}'))]); 
        [b,bint,resiMT{iSub},rint,stats] = regress(MTorg{iSub}', [PE_reg{iSub}' ones(size(PE_reg{iSub}'))]);         
    end
        
    % Preprocessing the decoded residuals (too large values should be truncated) 
    

%% Residual for uncertainties 
    load('UH_2Dim_voxelindex'); 
    nVox = length(ExtraFPvoxType); 
    for iSub = 1:21
        resiFP_byU{iSub} = zeros(size(DataMat_FP{iSub}.sfm)); 
        resiFP_byH{iSub} = zeros(size(DataMat_FP{iSub}.sfm)); 
        Uv = []; Hv = [];
        for iS = 1:4
            Uv = [Uv mDataMat.U_conv{iSub,iS}]; 
            Hv = [Hv mDataMat.H_conv{iSub,iS}]; 
        end
        for iVox = 1:nVox
            [~,~,resiFP_byU{iSub}(iVox,:),~,~] = regress(DataMat_FP{iSub}.sfm(iVox,:)', [Uv' ones(size(Uv'))]); 
            [~,~,resiFP_byH{iSub}(iVox,:),~,~] = regress(DataMat_FP{iSub}.sfm(iVox,:)', [Hv' ones(size(Hv'))]); 
        end 
    end
    
%     for iSub = 1:21
%         resiDec_re{iSub} = resiDec{iSub}; 
%         resiDec_re{iSub}(resiDec{iSub} > 2*std(resiDec{iSub})) = 2*std(resiDec{iSub}); 
%         resiDec_re{iSub}(resiDec{iSub} < -2*std(resiDec{iSub})) = -2*std(resiDec{iSub});         
%     end

    
    for iSub = 1:21
        % Calculate the correlation (explained variance) 
        cMt(iSub) = corr(resiDec_re{iSub}, resiMT{iSub}); 
        cMt_org(iSub) = corr(slopeDec{iSub}', MTv'); 
        cMt_cw(iSub,1) = corr(resiDec_re{iSub}(Dec_reg{iSub}==1), resiMT{iSub}(Dec_reg{iSub}==1)); 
        cMt_org_cw(iSub,1) = corr(slopeDec{iSub}(Dec_reg{iSub}==1)', MTv(Dec_reg{iSub}==1)'); 
        cMt_cw(iSub,2) = corr(resiDec_re{iSub}(Dec_reg{iSub}~=1), resiMT{iSub}(Dec_reg{iSub}~=1)); 
        cMt_org_cw(iSub,2) = corr(slopeDec{iSub}(Dec_reg{iSub}~=1)', MTv(Dec_reg{iSub}~=1)'); 
    end
    
    dt =  2.4;
    xRange = 0:dt:max(cHRF.t);
    for isamp = 1:length(xRange)
        sampleIndex = find(abs(xRange(isamp)-cHRF.t) == min(abs(xRange(isamp)-cHRF.t)),1,'first'); 
        cHRF.val_undersampled(isamp) = cHRF.val(sampleIndex); 
        cHRF.t_undersampled(isamp) = cHRF.t(sampleIndex); 
    end
    
    % Calculate the residual uncertainties
    for iSub = 1:21
%         PEtemp = abs(resiDec_re{iSub}-Percept_all{iSub}');
%         matEvent_conv = median(dt)*conv(cHRF.val_undersampled, PEtemp); 
%         extraU{iSub} = matEvent_conv(1:length(PEtemp))'; 
        extraU{iSub} = abs(resiDec{iSub}-Percept_all{iSub}')'; 
    end
    
    
    % Calculcate the ROI-wise residual correlation 
    for iROI = 1:2
        if iROI == 1
            valIndROI{iROI} = (DataMat.nClusterROI.assignments == 4) & (FPvoxType == 1); 
        else
            valIndROI{iROI} = (DataMat.nClusterROI.assignments == 13) & (FPvoxType == 2); 
        end
    end
    for iSub = 1:21
        for iROI = 1:2
            if iROI == 2
                resiFP_byUROI{iROI,iSub} = mean(resiFP_byU{iSub}(valIndROI{iROI},:)); 
            else
                resiFP_byUROI{iROI,iSub} = mean(resiFP_byH{iSub}(valIndROI{iROI},:)); 
            end
        end
    end
    
    % Calculate the residual correlations 
    for iSub = 1:21
        cFP_ROI(iSub) = corr(extraU{iSub}', resiFP_byUROI{2,iSub}');
        cFP_cwROI(1,iSub) = corr(extraU{iSub}(Dec_reg{iSub}==1)', resiFP_byUROI{2,iSub}(1,Dec_reg{iSub}==1)'); 
        cFP_cwROI(2,iSub) = corr(extraU{iSub}(Dec_reg{iSub}~=1)', resiFP_byUROI{2,iSub}(1,Dec_reg{iSub}~=1)'); 
%         for iVox = 1:nVox
%             cFP_cw{1}(iVox,iSub) = corr(extraU{iSub}(Dec_reg{iSub}==1)', resiFP_byU{iSub}(iVox,Dec_reg{iSub}==1)'); 
%             cFP_cw{2}(iVox,iSub) = corr(extraU{iSub}(Dec_reg{iSub}~=1)', resiFP_byU{iSub}(iVox,Dec_reg{iSub}~=1)'); 
%     %         cMt_org_cw(iSub,1) = corr(slopeDec{iSub}(Dec_reg{iSub}==1)', MTv(Dec_reg{iSub}==1)'); 
%     %         cFP_cw{2}(iVox,iSub) = corr(resiDec{iSub}(Dec_reg{iSub}~=1), resiFP_byU{iSub}(iVox,Dec_reg{iSub}~=1)'); 
%     %         cMt_org_cw(iSub,2) = corr(slopeDec{iSub}(Dec_reg{iSub}~=1)', MTv(Dec_reg{iSub}~=1)');
%         end
    end


%     [a1,a2] = ttest(cFP_cw{1}')
%     sum(a1==1)
%     [a1,a2] = ttest(cFP_cw{2}')
%     sum(a1==1)
    
    [a1,a2] = ttest(cFP_ROI)




    
    
    

figure(1); clf; 
for iSub = 1:21
    subplot(3,7,iSub); hold on; 
    plot(PE_reg{iSub}(Dec_reg{iSub}==1), slopeDec{iSub}(Dec_reg{iSub}==1),'r.'); 
    plot(PE_reg{iSub}(Dec_reg{iSub}==-1), slopeDec{iSub}(Dec_reg{iSub}==-1),'b.'); 
    corX(iSub) = corr(PE_reg{iSub}', slopeDec{iSub}'); 
    corXchoice(iSub,1) = corr(PE_reg{iSub}(Dec_reg{iSub}==1)', slopeDec{iSub}(Dec_reg{iSub}==1)'); 
    corXchoice(iSub,2) = corr(PE_reg{iSub}(Dec_reg{iSub}==-1)', slopeDec{iSub}(Dec_reg{iSub}==-1)');     
end

figure(2); clf; 
for iSub = 1:21
    subplot(3,7,iSub); hold on; 
    plot(resiDec{iSub}(Dec_reg{iSub}==1), resiMT{iSub}(Dec_reg{iSub}==1), 'r.');
    plot(resiDec{iSub}(Dec_reg{iSub}~=1), resiMT{iSub}(Dec_reg{iSub}~=1), 'b.');
end




Bin_size = 30;
Bin_dt = 0.02;
peall = []; 
for iSub = 1:21
    peall = [peall resiDec_re{iSub}'];
end
maxPE = max(abs(peall)); 
Rrange = -maxPE:Bin_dt:maxPE; 

clear binBOLD peBin nBins
for iSub = 1:21
    for iBin = 1:(length(Rrange)-Bin_size)
        indCW = ((Dec_reg{iSub} == 1) & (resiDec_re{iSub}' > Rrange(iBin)) & (resiDec_re{iSub}' <= Rrange(iBin+Bin_size))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (resiDec_re{iSub}' > Rrange(iBin)) & (resiDec_re{iSub}' <= Rrange(iBin+Bin_size)));
        peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
        if (sum(indCW)>1)
            binBOLD.cw(iSub, iBin) = mean(MTorg{iSub}(indCW)); 
        else
            binBOLD.cw(iSub, iBin) = nan; 
        end
        if (sum(indCCW)>1)
            binBOLD.ccw(iSub, iBin) = mean(MTorg{iSub}(indCCW)); 
        else
            binBOLD.ccw(iSub, iBin) = nan; 
        end
        binBOLD.all(iSub, iBin) = mean(MTorg{iSub}([indCW + indCCW]==1)); 
        nBins(iSub,iBin,1) = sum(indCW); 
        nBins(iSub,iBin,2) = sum(indCCW); 
    end
end

cwValInd = sum(isnan(binBOLD.cw))<1; 
ccwValInd = sum(isnan(binBOLD.ccw))<1; 

set(figure(3),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
pcw = polyfit(nanmean(peBin(:,cwValInd)), nanmean(binBOLD.cw(:,cwValInd)),1);
pccw = polyfit(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD.ccw(:,ccwValInd)),1);
% plot_patch_transparent(mean(peBin), mean(binBOLD.all), std(binBOLD.all,0,1),cmap_c(3,:),0.7); 
plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD.ccw(:,ccwValInd))-pccw(2), nanstd(binBOLD.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD.ccw(:,ccwValInd))-pccw(2)+nanstd(binBOLD.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD.ccw(:,ccwValInd))-pccw(2)-nanstd(binBOLD.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD.ccw(:,ccwValInd))-pccw(2), 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(binBOLD.cw(:,cwValInd))-pcw(2), nanstd(binBOLD.cw(:,cwValInd),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(binBOLD.cw(:,cwValInd))-pcw(2)+nanstd(binBOLD.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binBOLD.cw(:,cwValInd))-pcw(2)-nanstd(binBOLD.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binBOLD.cw(:,cwValInd))-pcw(2), 'color',cmap_c(3,:),'linewidth',1.5);
line([0 0],[-0.03 0.03],'linestyle','--','color','k'); 
xlim([-4 4]); ylim([-0.03 0.03]); 
xlabel('Expectation'); ylabel('Weighted BOLD'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-1.5 0 1.5], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')

set(figure(4),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD.ccw(:,ccwValInd)), nanstd(binBOLD.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD.ccw(:,ccwValInd))+nanstd(binBOLD.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD.ccw(:,ccwValInd))-nanstd(binBOLD.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD.ccw(:,ccwValInd)), 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(binBOLD.cw(:,cwValInd)), nanstd(binBOLD.cw(:,cwValInd),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(binBOLD.cw(:,cwValInd))+nanstd(binBOLD.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binBOLD.cw(:,cwValInd))-nanstd(binBOLD.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binBOLD.cw(:,cwValInd)), 'color',cmap_c(3,:),'linewidth',1.5);
line([0 0],[-0.03 0.03],'linestyle','--','color','k'); 
xlim([-3 3]); ylim([-0.04 0.04]); 
xlabel({'Extra Expectation','(Residual of weighted BOLD)'}); ylabel('Extra MT dv (LLR)'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.03 0 0.03], 'XTick', [-3 0 3], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')


set(figure(5),'position',[1 554 173 54]); clf; SP = subplot(1,1,1); hold on; 
plot(nanmean(peBin(:,cwValInd)), mean(nBins(:,cwValInd,1),1),'color',cmap_c(3,:),'linewidth',1.3);
plot(nanmean(peBin(:,ccwValInd)), mean(nBins(:,ccwValInd,2),1),'color',cmap_c(1,:),'linewidth',1.3);
xlim([-3 3]); ylim([0 70]);
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [0 70], 'XTick', [-3 0 3], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')






%% Binning for the FP activities 

% Calculate the residual's standard deviation 
for iSub = 1:21
    stdOMPFC(iSub) = std(resiDec_re{iSub}); 
end
for iROI = 1:2
    for iSub = 1:21
        stdv(iROI, iSub) = std(resiFP_byUROI{iROI,iSub}); 
    end
end

Bin_size = 30;
Bin_dt = 0.02;
peall = []; 
for iSub = 1:21
    peall = [peall resiDec_re{iSub}'];
end
maxPE = max(abs(peall)); 
Rrange = -maxPE:Bin_dt:maxPE; 

clear binFP peBin nBins
for iROI = 1:2
    for iSub = 1:21
        for iBin = 1:(length(Rrange)-Bin_size)
            indCW = ((Percept_all{iSub} == 1) & (resiDec_re{iSub}' > Rrange(iBin)) & (resiDec_re{iSub}' <= Rrange(iBin+Bin_size))); 
            indCCW = ((Percept_all{iSub} ~= 1) & (resiDec_re{iSub}' > Rrange(iBin)) & (resiDec_re{iSub}' <= Rrange(iBin+Bin_size)));
            peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
            if (sum(indCW)>1)
                binFP.cw(iROI, iSub, iBin) = mean(resiFP_byUROI{iROI,iSub}(indCW)); 
            else
                binFP.cw(iROI, iSub, iBin) = nan; 
            end
            if (sum(indCCW)>1)
                binFP.ccw(iROI, iSub, iBin) = mean(resiFP_byUROI{iROI,iSub}(indCCW)); 
            else
                binFP.ccw(iROI, iSub, iBin) = nan; 
            end
            binFP.all(iROI, iSub, iBin) = mean(resiFP_byUROI{iROI,iSub}([indCW + indCCW]==1)); 
            nBins(iSub,iBin,1) = sum(indCW); 
            nBins(iSub,iBin,2) = sum(indCCW); 
        end
    end
end

cwValInd = sum(isnan(squeeze(binFP.cw(1,:,:))))<1; 
ccwValInd = sum(isnan(squeeze(binFP.ccw(1,:,:))))<1; 


% Drawing a figure for IPS
set(figure(6),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(binFP.ccw(1, :,ccwValInd))), nanstd(squeeze(binFP.ccw(1, :,ccwValInd)),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(binFP.ccw(1,:,ccwValInd)))+nanstd(squeeze(binFP.ccw(1,:,ccwValInd)),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(binFP.ccw(1,:,ccwValInd)))-nanstd(squeeze(binFP.ccw(1,:,ccwValInd)),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(binFP.ccw(1,:,ccwValInd))), 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(squeeze(binFP.cw(1,:,cwValInd))), nanstd(squeeze(binFP.cw(1,:,cwValInd)),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(squeeze(binFP.cw(1,:,cwValInd)))+nanstd(squeeze(binFP.cw(1,:,cwValInd)),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(squeeze(binFP.cw(1,:,cwValInd)))-nanstd(squeeze(binFP.cw(1,:,cwValInd)),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(squeeze(binFP.cw(1,:,cwValInd))), 'color',cmap_c(3,:),'linewidth',1.5);
% line([0 0],[-0.03 0.03],'linestyle','--','color','k'); 
% xlim([-4 4]); ylim([-0.03 0.03]); 
% xlabel({'Extra Expectation','(Residual of weighted BOLD)'}); ylabel('MT dv (LLR)'); 
% set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.03 0 0.03], 'XTick', [-4 0 4], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')




set(figure(7),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(binFP.ccw(2, :,ccwValInd))), nanstd(squeeze(binFP.ccw(2, :,ccwValInd)),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(binFP.ccw(2,:,ccwValInd)))+nanstd(squeeze(binFP.ccw(2,:,ccwValInd)),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(binFP.ccw(2,:,ccwValInd)))-nanstd(squeeze(binFP.ccw(2,:,ccwValInd)),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(binFP.ccw(2,:,ccwValInd))), 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(squeeze(binFP.cw(2,:,cwValInd))), nanstd(squeeze(binFP.cw(2,:,cwValInd)),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(squeeze(binFP.cw(2,:,cwValInd)))+nanstd(squeeze(binFP.cw(2,:,cwValInd)),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(squeeze(binFP.cw(2,:,cwValInd)))-nanstd(squeeze(binFP.cw(2,:,cwValInd)),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(squeeze(binFP.cw(2,:,cwValInd))), 'color',cmap_c(3,:),'linewidth',1.5);


clear peBin 
peall = []; 
for iSub = 1:21
    peall = [peall (resiDec_re{iSub}*stdOMPFC(iSub))'];
end
maxPE = max(abs(peall)); 
Rrange = -maxPE:Bin_dt:maxPE; 
for iROI = 1:2
    for iSub = 1:21
        for iBin = 1:(length(Rrange)-Bin_size)
            indCW = ((Percept_all{iSub} == 1) & ((resiDec_re{iSub}*stdOMPFC(iSub))' > Rrange(iBin)) & ((resiDec_re{iSub}*stdOMPFC(iSub))' <= Rrange(iBin+Bin_size))); 
            indCCW = ((Percept_all{iSub} ~= 1) & ((resiDec_re{iSub}*stdOMPFC(iSub))' > Rrange(iBin)) & ((resiDec_re{iSub}*stdOMPFC(iSub))' <= Rrange(iBin+Bin_size)));
            peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
            if (sum(indCW)>1)
                fbinFP.cw(iROI, iSub, iBin) = mean(resiFP_byUROI{iROI,iSub}(indCW)); 
            else
                fbinFP.cw(iROI, iSub, iBin) = nan; 
            end
            if (sum(indCCW)>1)
                fbinFP.ccw(iROI, iSub, iBin) = mean(resiFP_byUROI{iROI,iSub}(indCCW)); 
            else
                fbinFP.ccw(iROI, iSub, iBin) = nan; 
            end
            fbinFP.all(iROI, iSub, iBin) = mean(resiFP_byUROI{iROI,iSub}([indCW + indCCW]==1)); 
            nBins(iSub,iBin,1) = sum(indCW); 
            nBins(iSub,iBin,2) = sum(indCCW); 
        end
    end
end


cwValInd = sum(isnan(squeeze(fbinFP.cw(1,:,:))))<2; 
ccwValInd = sum(isnan(squeeze(fbinFP.ccw(1,:,:))))<2; 

% Drawing a figure for IPS
set(figure(26),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(fbinFP.ccw(1, :,ccwValInd))), nanstd(squeeze(fbinFP.ccw(1, :,ccwValInd)),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(fbinFP.ccw(1,:,ccwValInd)))+nanstd(squeeze(fbinFP.ccw(1,:,ccwValInd)),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(fbinFP.ccw(1,:,ccwValInd)))-nanstd(squeeze(fbinFP.ccw(1,:,ccwValInd)),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(fbinFP.ccw(1,:,ccwValInd))), 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(squeeze(fbinFP.cw(1,:,cwValInd))), nanstd(squeeze(fbinFP.cw(1,:,cwValInd)),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(squeeze(fbinFP.cw(1,:,cwValInd)))+nanstd(squeeze(fbinFP.cw(1,:,cwValInd)),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(squeeze(fbinFP.cw(1,:,cwValInd)))-nanstd(squeeze(fbinFP.cw(1,:,cwValInd)),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(squeeze(fbinFP.cw(1,:,cwValInd))), 'color',cmap_c(3,:),'linewidth',1.5);
% line([0 0],[-0.03 0.03],'linestyle','--','color','k'); 
% xlim([-4 4]); ylim([-0.03 0.03]); 
% xlabel({'Extra Expectation','(Residual of weighted BOLD)'}); ylabel('MT dv (LLR)'); 
% set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.03 0 0.03], 'XTick', [-4 0 4], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')


set(figure(27),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(fbinFP.ccw(2, :,ccwValInd))), nanstd(squeeze(fbinFP.ccw(2, :,ccwValInd)),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(fbinFP.ccw(2,:,ccwValInd)))+nanstd(squeeze(fbinFP.ccw(2,:,ccwValInd)),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(fbinFP.ccw(2,:,ccwValInd)))-nanstd(squeeze(fbinFP.ccw(2,:,ccwValInd)),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(squeeze(fbinFP.ccw(2,:,ccwValInd))), 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(squeeze(fbinFP.cw(2,:,cwValInd))), nanstd(squeeze(fbinFP.cw(2,:,cwValInd)),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(squeeze(fbinFP.cw(2,:,cwValInd)))+nanstd(squeeze(fbinFP.cw(2,:,cwValInd)),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(squeeze(fbinFP.cw(2,:,cwValInd)))-nanstd(squeeze(fbinFP.cw(2,:,cwValInd)),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(squeeze(fbinFP.cw(2,:,cwValInd))), 'color',cmap_c(3,:),'linewidth',1.5);









































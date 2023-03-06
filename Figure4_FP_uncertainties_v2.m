clear all; 

addpath('/Users/joonwon/Desktop/Matlab_local/Brain_residualAnalysis_2018June/PleaseWork/NeatVersion_library'); 
addpath('/Users/joonwon/Desktop/Matlab_local/Brain_residualAnalysis_2018June/PleaseWork'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/Brain2Brain'); 


%% Load behavioral data 
Setting_param_bhv;

%% Generating the Model time-series (cHRF-convolved, raw)
Generating_model_TS; 

% load('TS_OMPFC_all.mat');
% load('TS_MT_all.mat');
load('TS_FP_all15.mat');


PE_reg = cell(1,21); Dec_reg = cell(1,21); 
for iSub = 1:21
    for iS = 1:4
        PE_reg{iSub} = [PE_reg{iSub} mDataMat.PE_conv{iSub,iS}]; 
        Dec_reg{iSub} = [Dec_reg{iSub} mDataMat.Decision_conv{iSub,iS}];         
    end
end


%% Normalizing each bin's time-series
for iSub = 1:21
    DataMat_FP{iSub}.sfm_norm = (DataMat_FP{iSub}.sfm - repmat(mean(DataMat_FP{iSub}.sfm,2), 1, 500))./repmat(std(DataMat_FP{iSub}.sfm,0,2), 1,500);
end



%% 
load('ROI_coord_merge_25.mat','coordMat');







Bin_size = 25;
Bin_dt = 0.001;
peall = []; 
for iSub = 1:21
    peall = [peall PE_reg{iSub}];
end
maxPE = max(abs(peall)); 
Rrange = -maxPE:Bin_dt:maxPE; 




nBin_rd = 8; 
Rrange = linspace(-maxPE, maxPE, nBin_rd); 


clear binFP peBin nBins
binFP.cw = nan(21, length(DataMat_FP{1}.sfm_norm(:,1)), nBin_rd-1); 
binFP.ccw = nan(21, length(DataMat_FP{1}.sfm_norm(:,1)), nBin_rd-1); 
for iSub = 1:21
    iSub
    for iBin = 1:(nBin_rd-1)
        indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+1))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+1)));
        peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+1))/2; 
        for iVox = 1:length(DataMat_FP{iSub}.sfm_norm(:,1))
            if (sum(indCW)>1)
                binFP.cw(iSub, iVox, iBin) = mean(DataMat_FP{iSub}.sfm_norm(iVox,indCW)); 
            end
            if (sum(indCCW)>1)
                binFP.ccw(iSub, iVox, iBin) = mean(DataMat_FP{iSub}.sfm_norm(iVox,indCCW));
            end
            binFP.all(iSub, iVox, iBin) = mean(DataMat_FP{iSub}.sfm_norm(iVox,[indCW + indCCW]==1)); 
        end
    end
end

cwValInd = sum(isnan(squeeze(binFP.cw(:,1,:))))<2; 
ccwValInd = sum(isnan(squeeze(binFP.ccw(:,1,:))))<2; 


% Substitute nan value
for iSub = 1:21
    indVal = ~isnan(squeeze(binFP.ccw(iSub,1,:))); 
    for iVox = 1:length(binFP.ccw(1,:,1))
        pt = polyfit(nanmean(peBin(:,(indVal))), squeeze(binFP.ccw(iSub, iVox, (indVal)))', 1); 
        binFP.ccw(iSub, iVox, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
    end
    indVal = ~isnan(squeeze(binFP.cw(iSub,1,:)));
    for iVox = 1:length(binFP.cw(1,:,1))
        pt = polyfit(nanmean(peBin(:,(indVal))), squeeze(binFP.cw(iSub, iVox, (indVal)))', 1); 
        binFP.cw(iSub, iVox, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
    end
    indVal = ~isnan(squeeze(binFP.all(iSub,1,:)));
    for iVox = 1:length(binFP.all(1,:,1))
        pt = polyfit(nanmean(peBin(:,(indVal))), squeeze(binFP.all(iSub, iVox, (indVal)))', 1); 
        binFP.all(iSub, iVox, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
    end
end



TSa = squeeze(nanmean(binFP.all,1)); 
TSa_cw = squeeze(nanmean(binFP.cw,1)); 
TSa_ccw = squeeze(nanmean(binFP.ccw,1)); 

xind = nanmean(peBin);


% Discriminate the H- and U-preferring voxels 
% First condition: U vs. H (shape difference)
for iVox = 1:length(binFP.cw(1,:,1))
    [~,~,r1] = regress(TSa(iVox,:)', [(xind)' ones(size(xind))']); 
    [~,~,r2] = regress(TSa(iVox,:)', [(xind.^2)' ones(size(xind))']); 
    rmse.U_linear(iVox) = sqrt(sum(r1.^2)/length(peBin(1,:))); 
    rmse.H_symm(iVox) = sqrt(sum(r2.^2)/length(peBin(1,:))); 
end
ShapeDiff = - rmse.H_symm + rmse.U_linear ; 

% Second condition: U vs. H (choice difference)
for iVox = 1:length(binFP.cw(1,:,1))
    ChoiceDiff(iVox) = corr(TSa_cw(iVox,:)', TSa_ccw(iVox,:)'); 
end

figure(61); clf; 
SP = subplot(1,1,1); cla; hold on; 
plot(ChoiceDiff(ChoiceDiff>0 & ShapeDiff>0), ShapeDiff(ChoiceDiff>0 & ShapeDiff>0),'k.','MarkerFaceColor','k','MarkerEdgeColor','k','marker','.'); 
plot(ChoiceDiff(ChoiceDiff<0 & ShapeDiff<0), ShapeDiff(ChoiceDiff<0 & ShapeDiff<0),'k.','MarkerFaceColor','k','MarkerEdgeColor','k','marker','.'); 
plot(ChoiceDiff(ChoiceDiff<0 & ShapeDiff>0), ShapeDiff(ChoiceDiff<0 & ShapeDiff>0),'k.','MarkerFaceColor',[0 0 0]+0.7,'MarkerEdgeColor',[0 0 0]+0.7,'marker','.'); 
plot(ChoiceDiff(ChoiceDiff>0 & ShapeDiff<0), ShapeDiff(ChoiceDiff>0 & ShapeDiff<0),'k.','MarkerFaceColor',[0 0 0]+0.7,'MarkerEdgeColor',[0 0 0]+0.7,'marker','.'); 
line([-1 1],[0 0],'linestyle','--','color','k'); 
line([0 0],[-0.2 0.2],'linestyle','--','color','k'); 
xlim([-1 1]); ylim([-0.2 0.2]); 


xlabel({'Between-choice', 'pattern similarity'}); ylabel('Shape of function'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.2 0 0.2], 'XTick', [-1 0 1], 'FontSize', 10,'color','none')



nVox = length(binFP.cw(1,:,1));

nBin = 8; 
clear BinFP BinPE
BinPE = nan(21,nBin-1); 
BinFP = nan(21,nVox,nBin-1); 
Rrange = linspace(-maxPE, maxPE, nBin_rd); 

for iSub = 1:21
    % Divide the time-range
%     Rrange = linspace(min(PE_reg{iSub}), max(PE_reg{iSub}),nBin); 
    BinPE(iSub,:) = nan(1,nBin-1); 
    

    for iBin = 1:(nBin-1)
        peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+1))/2; 
        indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+1))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(nBin-iBin)) & (PE_reg{iSub} <= Rrange(nBin-iBin+1))); 
    
        BinPE(iSub,iBin) = nanmean([PE_reg{iSub}(indCW) -PE_reg{iSub}(indCCW)]); 
        
        for iVox = 1:nVox
            BinFP(iSub,iVox,iBin) = nanmean([DataMat_FP{iSub}.sfm(iVox,indCW) DataMat_FP{iSub}.sfm(iVox,indCCW)]); 
        end
    end
end
TSa = squeeze(nanmean(BinFP,1)); 
clear fpBinFP fpBinPE
fpBinPE{1} = nan(21,nBin-1); 
fpBinPE{2} = nan(21,nBin-1); 
fpBinFP{1} = nan(21,nVox,nBin-1); 
fpBinFP{2} = nan(21,nVox,nBin-1); 

for iSub = 1:21
    % Divide the time-range
%     Rrange = linspace(min(PE_reg{iSub}), max(PE_reg{iSub}),nBin); 
    fpBinPE{1}(iSub,:) = nan(1,nBin-1); 
    fpBinPE{2}(iSub,:) = nan(1,nBin-1); 
    for iBin = 1:(nBin-1)
        indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+1))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+1))); 
    
        fpBinPE{1}(iSub,iBin) = nanmean(PE_reg{iSub}(indCW)); 
        fpBinPE{2}(iSub,iBin) = nanmean(PE_reg{iSub}(indCCW)); 
        
        for iVox = 1:nVox
            if (sum(indCW)>1)
                fpBinFP{1}(iSub,iVox,iBin) = nanmean(DataMat_FP{iSub}.sfm(iVox,indCW)); 
            else
                fpBinFP{1}(iSub,iVox,iBin) = nan; 
            end
            if (sum(indCCW)>1)
                fpBinFP{2}(iSub,iVox,iBin) = nanmean(DataMat_FP{iSub}.sfm(iVox,indCCW)); 
            else
                fpBinFP{2}(iSub,iVox,iBin) = nan;
            end
        end
    end
end
% Substitute nan value
for iSub = 1:21
    indVal = ~isnan(squeeze(fpBinFP{1}(iSub,1,:))); 
    for iVox = 1:length(binFP.ccw(1,:,1))
        pt = polyfit(nanmean(peBin(:,(indVal))), squeeze(fpBinFP{1}(iSub, iVox, (indVal)))', 1); 
        fpBinFP{1}(iSub, iVox, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
    end
    
    indVal = ~isnan(squeeze(fpBinFP{2}(iSub,1,:)));
    for iVox = 1:length(binFP.cw(1,:,1))
        pt = polyfit(nanmean(peBin(:,(indVal))), squeeze(fpBinFP{2}(iSub, iVox, (indVal)))', 1); 
        fpBinFP{2}(iSub, iVox, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
    end
end



TSa_all = TSa;
TSa_cw = squeeze(nanmean(fpBinFP{1},1)); 
TSa_ccw = squeeze(nanmean(fpBinFP{2},1)); 

xind = nanmean(BinPE);

% First condition: U vs. H (shape difference)
for iVox = 1:nVox
    [~,~,r1] = regress(TSa(iVox,:)', [(xind)' ones(size(xind))']); 
    [~,~,r2] = regress(TSa(iVox,:)', [(xind.^2)' ones(size(xind))']); 
    rmse.U_linear(iVox) = sqrt(sum(r1.^2)/(nBin-2)); 
    rmse.H_symm(iVox) = sqrt(sum(r2.^2)/(nBin-2)); 
end
ShapeDiff = - rmse.H_symm + rmse.U_linear ; 

% Second condition: U vs. H (choice difference)
for iVox = 1:nVox
    ChoiceDiff(iVox) = corr(TSa_cw(iVox,2:7)', TSa_ccw(iVox,:)'); 
end

clear FPvoxType ExtraFPvoxType; 
for iVox = 1:nVox
    if (ChoiceDiff(iVox) > 0) && (ShapeDiff(iVox) > 0)
        FPvoxType(iVox) = 1; % H-preferring
    elseif (ChoiceDiff(iVox) < 0) && (ShapeDiff(iVox) < 0)
        FPvoxType(iVox) = 2; % U-preferring
    else
        FPvoxType(iVox) = 3; % Else
    end
end
for iVox = 1:nVox
    if (ChoiceDiff(iVox) > 0.5) && (ShapeDiff(iVox) > 0.05)
        ExtraFPvoxType(iVox) = 1; % H-preferring
    elseif (ChoiceDiff(iVox) < -0.5) && (ShapeDiff(iVox) < -0.050)
        ExtraFPvoxType(iVox) = 2; % U-preferring
    else
        ExtraFPvoxType(iVox) = 3; % Else
    end
end

figure(62); clf; 
SP = subplot(1,1,1); cla; hold on; 
plot(ChoiceDiff(ChoiceDiff>0 & ShapeDiff>0), ShapeDiff(ChoiceDiff>0 & ShapeDiff>0),'k.','MarkerFaceColor','k','MarkerEdgeColor','k','marker','.'); 
plot(ChoiceDiff(ChoiceDiff<0 & ShapeDiff<0), ShapeDiff(ChoiceDiff<0 & ShapeDiff<0),'k.','MarkerFaceColor','k','MarkerEdgeColor','k','marker','.'); 
plot(ChoiceDiff(ChoiceDiff<0 & ShapeDiff>0), ShapeDiff(ChoiceDiff<0 & ShapeDiff>0),'k.','MarkerFaceColor',[0 0 0]+0.7,'MarkerEdgeColor',[0 0 0]+0.7,'marker','.'); 
plot(ChoiceDiff(ChoiceDiff>0 & ShapeDiff<0), ShapeDiff(ChoiceDiff>0 & ShapeDiff<0),'k.','MarkerFaceColor',[0 0 0]+0.7,'MarkerEdgeColor',[0 0 0]+0.7,'marker','.'); 
line([-1 1],[0 0],'linestyle','--','color','k'); 
line([0 0],[-0.2 0.2],'linestyle','--','color','k'); 
xlim([-1 1]); ylim([-0.2 0.2]); 


xlabel({'Between-choice', 'pattern similarity'}); ylabel('Shape of function'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.2 0 0.2], 'XTick', [-1 0 1], 'FontSize', 10,'color','none')




save('UH_2Dim_voxelindex.mat','FPvoxType','ExtraFPvoxType'); 
load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7/Data_allTS_15.mat'); 

for iROI = 1:2
    if iROI == 1
        valIndROI{iROI} = (DataMat.nClusterROI.assignments == 4) & (FPvoxType == 1); 
    else
        valIndROI{iROI} = (DataMat.nClusterROI.assignments == 13) & (FPvoxType == 2); 
    end
end


% 
clear fbinFP peBin nBins
Bin_size = 20;
Bin_dt = 0.001;
peall = []; 
for iSub = 1:21
    peall = [peall PE_reg{iSub}];
end
maxPE = max(abs(peall)); 
Rrange = -maxPE:Bin_dt:maxPE; 
for iROI = 1:2
    fbinFP{iROI}.cw = nan(21, nBin_rd-1); 
    fbinFP{iROI}.ccw = nan(21, nBin_rd-1); 
    fbinFP{iROI}.all = nan(21, nBin_rd-1); 
    for iSub = 1:21
        for iBin = 1:(nBin_rd-1)
            indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size))); 
            indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size)));
            peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
            
            fbinFP{iROI}.cw(iSub,iBin) = mean(mean(DataMat_FP{iSub}.sfm_norm(valIndROI{iROI},indCW))); 
            fbinFP{iROI}.ccw(iSub,iBin) = mean(mean(DataMat_FP{iSub}.sfm_norm(valIndROI{iROI},indCCW))); 
            fbinFP{iROI}.all(iSub,iBin) = mean(mean(DataMat_FP{iSub}.sfm_norm(valIndROI{iROI},[indCW+indCCW]==1))); 
        end
    end
end

cwValInd = sum(isnan(squeeze(fbinFP{1}.cw)))<2; 
ccwValInd = sum(isnan(squeeze(fbinFP{1}.ccw)))<2; 


% % Substitute nan value
% for iROI = 1:2
%     for iSub = 1:21
%         indVal = ~isnan(squeeze(fbinFP{iROI}.ccw(iSub,:))); 
%         pt = polyfit(nanmean(peBin(:,(indVal))), squeeze(fbinFP{iROI}.ccw(iSub, (indVal))), 1); 
%         fbinFP{iROI}.ccw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
%         
%         indVal = ~isnan(squeeze(fbinFP{iROI}.cw(iSub,:))); 
%         pt = polyfit(nanmean(peBin(:,(indVal))), squeeze(fbinFP{iROI}.cw(iSub, (indVal))), 1); 
%         fbinFP{iROI}.cw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
%         
%         indVal = ~isnan(squeeze(fbinFP{iROI}.all(iSub,:))); 
%         pt = polyfit(nanmean(peBin(:,(indVal))), squeeze(fbinFP{iROI}.all(iSub, (indVal))), 1); 
%         fbinFP{iROI}.all(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
%     end
% end

set(figure(208),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
iROI = 1; 
plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP{iROI}.ccw(:,ccwValInd)), nanstd(fbinFP{iROI}.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP{iROI}.ccw(:,ccwValInd))+nanstd(fbinFP{iROI}.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP{iROI}.ccw(:,ccwValInd))-nanstd(fbinFP{iROI}.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP{iROI}.ccw(:,ccwValInd)), 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(fbinFP{iROI}.cw(:,cwValInd)), nanstd(fbinFP{iROI}.cw(:,cwValInd),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(fbinFP{iROI}.cw(:,cwValInd))+nanstd(fbinFP{iROI}.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(fbinFP{iROI}.cw(:,cwValInd))-nanstd(fbinFP{iROI}.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(fbinFP{iROI}.cw(:,cwValInd)), 'color',cmap_c(3,:),'linewidth',1.5);
line([0 0],[-0.2 0.6],'linestyle','--','color','k'); 
xlim([-0.05 0.05]); ylim([-0.2 0.6]); 
xlabel('Expectation'); ylabel('IPS BOLD'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.2:0.2:0.6], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')


set(figure(209),'position',[1 212 173 161]); clf; SP = subplot(1,1,1); hold on; 
iROI = 2; 
plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP{iROI}.ccw(:,ccwValInd)), nanstd(fbinFP{iROI}.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP{iROI}.ccw(:,ccwValInd))+nanstd(fbinFP{iROI}.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP{iROI}.ccw(:,ccwValInd))-nanstd(fbinFP{iROI}.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP{iROI}.ccw(:,ccwValInd)), 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(fbinFP{iROI}.cw(:,cwValInd)), nanstd(fbinFP{iROI}.cw(:,cwValInd),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(fbinFP{iROI}.cw(:,cwValInd))+nanstd(fbinFP{iROI}.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(fbinFP{iROI}.cw(:,cwValInd))-nanstd(fbinFP{iROI}.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(fbinFP{iROI}.cw(:,cwValInd)), 'color',cmap_c(3,:),'linewidth',1.5);
line([0 0],[-0.2 0.6],'linestyle','--','color','k'); 
xlim([-0.05 0.05]); ylim([-0.2 0.6]); 
xlabel('Expectation'); ylabel('IPS BOLD'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.2:0.2:0.6], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')







set(figure(209),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
plot(nanmean(peBin(:,cwValInd)), mean(nBins(:,cwValInd,1),1),'color',cmap_c(3,:),'linewidth',1.3);
plot(nanmean(peBin(:,ccwValInd)), mean(nBins(:,ccwValInd,2),1),'color',cmap_c(1,:),'linewidth',1.3);
xlim([-0.05 0.05]); ylim([0 200]);
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [0 200], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')

















































































































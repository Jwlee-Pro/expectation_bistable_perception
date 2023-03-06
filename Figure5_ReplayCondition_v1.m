

clear all; 

addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2018_Mfin/Brain_residualAnalysis_2018June/PleaseWork/NeatVersion_library'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/Brain2Brain'); 


%% Load behavioral data 
Setting_param_bhv;

%% Generating the Model time-series (cHRF-convolved, raw)
Generating_model_TS; 

load('TS_OMPFC_all.mat');
load('TS_MT_all.mat');
load('TS_FP_all15.mat');

nSess = [ones(1,14) 4*ones(1,7)]; 


PE_reg = cell(1,21); Dec_reg = cell(1,21); 
PE_regorg = cell(1,21); Dec_regorg = cell(1,21); 
for iSub = 1:21
    for iS = 1:nSess(iSub)
        PE_reg{iSub} = [PE_reg{iSub} rDataMat.PE_conv{iSub,iS}]; 
        Dec_reg{iSub} = [Dec_reg{iSub} rDataMat.Decision_conv{iSub,iS}];         
    end
    for iS = 1:4
        PE_regorg{iSub} = [PE_regorg{iSub} mDataMat.PE_conv{iSub,iS}]; 
        Dec_regorg{iSub} = [Dec_regorg{iSub} mDataMat.Decision_conv{iSub,iS}];         
    end
    CorrVal.OMPFC_PEmimic{iSub} = corr(DataMat_OMPFC{iSub}.sfm', PE_regorg{iSub}');
    CorrVal.MT_Dmimic{iSub} = corr(DataMat_MT{iSub}.sfm', Dec_regorg{iSub}');
    CorrVal.OMPFC_PE{iSub} = corr(DataMat_OMPFC{iSub}.mimic', PE_reg{iSub}');
    CorrVal.MT_D{iSub} = corr(DataMat_MT{iSub}.mimic', Dec_reg{iSub}');
end

for iSub = 1:21
    DataMat_FP{iSub}.mimic_norm = (DataMat_FP{iSub}.mimic - repmat(mean(DataMat_FP{iSub}.mimic,2), 1, length(PE_reg{iSub})))./repmat(std(DataMat_FP{iSub}.mimic,0,2), 1,length(PE_reg{iSub}));
    DataMat_OMPFC{iSub}.mimic_norm = (DataMat_OMPFC{iSub}.mimic - repmat(mean(DataMat_OMPFC{iSub}.mimic,2), 1, length(PE_reg{iSub})))./repmat(std(DataMat_OMPFC{iSub}.mimic,0,2), 1,length(PE_reg{iSub}));
end


nBin = 6; 
clear BinFP BinPE
BinPE = nan(21,nBin-1); 
nVox = length(DataMat_FP{1}.sfm(:,1)); 
BinFP = nan(21,nVox,nBin-1); 

for iSub = 1:21
    % Divide the time-range
    Rrange = linspace(min(PE_reg{iSub}), max(PE_reg{iSub}),nBin); 
    BinPE(iSub,:) = nan(1,nBin-1); 
    for iBin = 1:(nBin-1)
        indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+1))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(nBin-iBin)) & (PE_reg{iSub} <= Rrange(nBin-iBin+1))); 
    
        BinPE(iSub,iBin) = nanmean([PE_reg{iSub}(indCW) -PE_reg{iSub}(indCCW)]); 
        
        for iVox = 1:nVox
            BinFP(iSub,iVox,iBin) = nanmean([DataMat_FP{iSub}.mimic_norm(iVox,indCW) DataMat_FP{iSub}.mimic_norm(iVox,indCCW)]); 
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
    Rrange = linspace(min(PE_reg{iSub}), max(PE_reg{iSub}),nBin); 
    fpBinPE{1}(iSub,:) = nan(1,nBin-1); 
    fpBinPE{2}(iSub,:) = nan(1,nBin-1); 
    for iBin = 1:(nBin-1)
        indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+1))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+1))); 
    
        fpBinPE{1}(iSub,iBin) = nanmean(PE_reg{iSub}(indCW)); 
        fpBinPE{2}(iSub,iBin) = nanmean(PE_reg{iSub}(indCCW)); 
        
        for iVox = 1:nVox
            fpBinFP{1}(iSub,iVox,iBin) = nanmean(DataMat_FP{iSub}.mimic_norm(iVox,indCW)); 
            fpBinFP{2}(iSub,iVox,iBin) = nanmean(DataMat_FP{iSub}.mimic_norm(iVox,indCCW)); 
        end
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
    ChoiceDiff(iVox) = corr(TSa_cw(iVox,:)', TSa_ccw(iVox,:)'); 
end


for iVox = 1:nVox
    if (ChoiceDiff(iVox) > 0) && (ShapeDiff(iVox) > 0)
        FPvoxType_mimic(iVox) = 1; % H-preferring
    elseif (ChoiceDiff(iVox) < 0) && (ShapeDiff(iVox) < 0)
        FPvoxType_mimic(iVox) = 2; % U-preferring
    else
        FPvoxType_mimic(iVox) = 3; % Else
    end
end

for iVox = 1:nVox
    if (ChoiceDiff(iVox) > 0.5) && (ShapeDiff(iVox) > 0.025)
        ExtraFPvoxType_mimic(iVox) = 1; % H-preferring
    elseif (ChoiceDiff(iVox) < -0.5) && (ShapeDiff(iVox) < -0.0250)
        ExtraFPvoxType_mimic(iVox) = 2; % U-preferring
    else
        ExtraFPvoxType_mimic(iVox) = 3; % Else
    end
end








load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7/Data_allTS_15.mat'); 
load('UH_2Dim_voxelindex'); 


for iROI = 1:2
    if iROI == 1
        valIndROI_re{iROI} = (DataMat.nClusterROI.assignments == 4) & (ExtraFPvoxType_mimic == 1); 
    else
        valIndROI_re{iROI} = (DataMat.nClusterROI.assignments == 7) & (ExtraFPvoxType_mimic == 2); 
    end
end



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
    fbinFP{iROI}.cw = nan(21, (length(Rrange)-Bin_size)); 
    fbinFP{iROI}.ccw = nan(21, (length(Rrange)-Bin_size)); 
    fbinFP{iROI}.all = nan(21, (length(Rrange)-Bin_size)); 
    for iSub = 1:21
        for iBin = 1:(length(Rrange)-Bin_size)
            indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size))); 
            indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size)));
            peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
            
            fbinFP{iROI}.cw(iSub,iBin) = mean(mean(DataMat_FP{iSub}.mimic_norm(valIndROI_re{iROI},indCW))); 
            fbinFP{iROI}.ccw(iSub,iBin) = mean(mean(DataMat_FP{iSub}.mimic_norm(valIndROI_re{iROI},indCCW))); 
            fbinFP{iROI}.all(iSub,iBin) = mean(mean(DataMat_FP{iSub}.mimic_norm(valIndROI_re{iROI},[indCW+indCCW]==1))); 
        end
    end
end

cwValInd = sum(isnan(squeeze(fbinFP{1}.cw)))<6; 
ccwValInd = sum(isnan(squeeze(fbinFP{1}.ccw)))<6; 


% Substitute nan value
for iROI = 1:2
    for iSub = 1:21
        indVal = ~isnan(squeeze(fbinFP{iROI}.ccw(iSub,:))); 
        pt = polyfit(nanmean(peBin(:,(indVal))), squeeze(fbinFP{iROI}.ccw(iSub, (indVal))), 1); 
        fbinFP{iROI}.ccw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
        
        indVal = ~isnan(squeeze(fbinFP{iROI}.cw(iSub,:))); 
        pt = polyfit(nanmean(peBin(:,(indVal))), squeeze(fbinFP{iROI}.cw(iSub, (indVal))), 1); 
        fbinFP{iROI}.cw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
        
        indVal = ~isnan(squeeze(fbinFP{iROI}.all(iSub,:))); 
        pt = polyfit(nanmean(peBin(:,(indVal))), squeeze(fbinFP{iROI}.all(iSub, (indVal))), 1); 
        fbinFP{iROI}.all(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
    end
end













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

plot(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP_main{iROI}.ccw(:,ccwValInd)),'k-', 'color',cmap_choiceStrong(1,:),'linewidth',2);
plot(nanmean(peBin(:,cwValInd)), nanmean(fbinFP_main{iROI}.cw(:,cwValInd)),'k-', 'color',cmap_choiceStrong(3,:),'linewidth',2);

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
plot(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP_main{iROI}.ccw(:,ccwValInd)),'k-', 'color',cmap_choiceStrong(1,:),'linewidth',2);
plot(nanmean(peBin(:,cwValInd)), nanmean(fbinFP_main{iROI}.cw(:,cwValInd)),'k-', 'color',cmap_choiceStrong(3,:),'linewidth',2);

line([0 0],[-0.2 0.6],'linestyle','--','color','k'); 
xlim([-0.05 0.05]); ylim([-0.2 0.6]); 
xlabel('Expectation'); ylabel('IPS BOLD'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.2:0.2:0.6], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')




% SfM condition 

PE_reg = cell(1,21); Dec_reg = cell(1,21); 
for iSub = 1:21
    for iS = 1:4
        PE_reg{iSub} = [PE_reg{iSub} mDataMat.PE_conv{iSub,iS}]; 
        Dec_reg{iSub} = [Dec_reg{iSub} mDataMat.Decision_conv{iSub,iS}];         
    end
end
for iSub = 1:21
    DataMat_FP{iSub}.sfm_norm = (DataMat_FP{iSub}.sfm - repmat(mean(DataMat_FP{iSub}.sfm,2), 1, 500))./repmat(std(DataMat_FP{iSub}.sfm,0,2), 1,500);
end

for iROI = 1:2
    if iROI == 1
        valIndROI{iROI} = (DataMat.nClusterROI.assignments == 4) & (FPvoxType == 1); 
    else
        valIndROI{iROI} = (DataMat.nClusterROI.assignments == 13) & (FPvoxType == 2); 
    end
end


nBin_rd = 8; 
% Rrange = linspace(-maxPE, maxPE, nBin_rd); 

maxPE = max(abs(peall)); 
Rrange = -maxPE:Bin_dt:maxPE; 
for iROI = 1:2
    fbinFP_main{iROI}.cw = nan(21, (length(Rrange)-Bin_size)); 
    fbinFP_main{iROI}.ccw = nan(21, (length(Rrange)-Bin_size)); 
    fbinFP_main{iROI}.all = nan(21, (length(Rrange)-Bin_size)); 
    for iSub = 1:21
        for iBin = 1:(length(Rrange)-Bin_size)%(nBin_rd-1)
            indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size))); 
            indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size)));
            peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
            
            fbinFP_main{iROI}.cw(iSub,iBin) = mean(mean(DataMat_FP{iSub}.sfm_norm(valIndROI{iROI},indCW))); 
            fbinFP_main{iROI}.ccw(iSub,iBin) = mean(mean(DataMat_FP{iSub}.sfm_norm(valIndROI{iROI},indCCW))); 
            fbinFP_main{iROI}.all(iSub,iBin) = mean(mean(DataMat_FP{iSub}.sfm_norm(valIndROI{iROI},[indCW+indCCW]==1))); 
        end
    end
end

cwValInd = sum(isnan(squeeze(fbinFP_main{1}.cw)))<2; 
ccwValInd = sum(isnan(squeeze(fbinFP_main{1}.ccw)))<2; 



set(figure(308),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
iROI = 1; 
plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP_main{iROI}.ccw(:,ccwValInd)), nanstd(fbinFP_main{iROI}.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP_main{iROI}.ccw(:,ccwValInd))+nanstd(fbinFP_main{iROI}.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP_main{iROI}.ccw(:,ccwValInd))-nanstd(fbinFP_main{iROI}.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(fbinFP_main{iROI}.ccw(:,ccwValInd)), 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(fbinFP_main{iROI}.cw(:,cwValInd)), nanstd(fbinFP_main{iROI}.cw(:,cwValInd),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(fbinFP_main{iROI}.cw(:,cwValInd))+nanstd(fbinFP_main{iROI}.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(fbinFP_main{iROI}.cw(:,cwValInd))-nanstd(fbinFP_main{iROI}.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(fbinFP_main{iROI}.cw(:,cwValInd)), 'color',cmap_c(3,:),'linewidth',1.5);
line([0 0],[-0.2 0.6],'linestyle','--','color','k'); 
xlim([-0.05 0.05]); ylim([-0.2 0.6]); 
xlabel('Expectation'); ylabel('IPS BOLD'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.2:0.2:0.6], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')



figure(38); clf; hold on; 
plot(nanmean(fbinFP_main{iROI}.cw), nanmean(fbinFP{iROI}.cw), 'ko','markerfacecolor',cmap_c(3,:),'markeredgecolor','none','markersize',4); 
plot(nanmean(fbinFP_main{iROI}.ccw), nanmean(fbinFP{iROI}.ccw), 'ko','markerfacecolor',cmap_c(1,:),'markeredgecolor','none','markersize',4); 
line([-1 1],[-1 1],'linestyle','--','color','k');
% line([-1 1],[1 -1],'linestyle','--','color','k');
xlim([-1 1]); ylim([-1 1]);


% a = [nanmean(fbinFP_main{1}.cw) nanmean(fbinFP_main{1}.ccw)]; 
a = [nanmean(fbinFP_main{1}.cw)-nanmean(fbinFP_main{1}.cw(:,39)) nanmean(fbinFP_main{1}.ccw)-nanmean(fbinFP_main{1}.ccw(:,39))]; 
% b = [nanmean(fbinFP{1}.ccw) nanmean(fbinFP{1}.ccw)]; 
b = [nanmean(fbinFP{1}.cw)-nanmean(fbinFP{1}.cw(:,39)) nanmean(fbinFP{1}.ccw)-nanmean(fbinFP{1}.ccw(:,39))]; 
% b = b - nanmean(b);
[h,p,ci,stats] = ttest(a-b)

a = [nanmean(fbinFP_main{2}.cw) nanmean(fbinFP_main{2}.ccw)]; 
% a = a - nanmean(a); 
b = [nanmean(fbinFP{2}.ccw) nanmean(fbinFP{2}.ccw)]; 
% b = b - nanmean(b);
[h,p,ci,stats] = ttest(a-b)


%% MT 

% Load decoded decision variable from MT 
% MT data decoding 
load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7/MTanalysis_deep/LeeLee_PriorEffect_workspace_PredCo_2018June.mat'); 

MTdec = cell(1,21); 
MTdec_mimic = cell(1,21); 
for iSub = 1:21
    for iS = 1:4
        ta = dec_final_fin{5}{iSub,iS}; 
%         ta = dec_final_fin{4}{iSub,iS}; 
        MTdec{iSub} = [MTdec{iSub} log(ta)-log(1-ta)];
    end
    for iS = 1:nSess(iSub)
%         ta = dec_final_fin_mimic{4}{iSub,iS}; 
        ta = dec_final_fin_mimic{5}{iSub,iS}; 
        MTdec_mimic{iSub} = [MTdec_mimic{iSub} log(ta)-log(1-ta)]; 
    end
end

clear binMTdv
for iSub = 1:21
    for iBin = 1:(length(Rrange)-Bin_size)
        indCW = ((Dec_regorg{iSub} == 1) & (PE_regorg{iSub} > Rrange(iBin)) & (PE_regorg{iSub} <= Rrange(iBin+Bin_size))); 
        indCCW = ((Dec_regorg{iSub} ~= 1) & (PE_regorg{iSub} > Rrange(iBin)) & (PE_regorg{iSub} <= Rrange(iBin+Bin_size)));
        if (sum(indCW)>=1)
            binMTdv.cw(iSub, iBin) = mean(MTdec{iSub}(indCW)); 
        else
            binMTdv.cw(iSub, iBin) = nan; 
        end
        if (sum(indCCW)>=1)
            binMTdv.ccw(iSub, iBin) = mean(MTdec{iSub}(indCCW)); 
        else
            binMTdv.ccw(iSub, iBin) = nan; 
        end
        binMTdv.all(iSub, iBin) = mean(MTdec{iSub}([indCW + indCCW]==1)); 
    end
end

cwValInd = sum(isnan(binMTdv.cw))<2; 
ccwValInd = sum(isnan(binMTdv.ccw))<2; 


% Substitute nan value
for iSub = 1:21
    indVal = ~isnan(binMTdv.ccw(iSub,:)); 
    pt = polyfit(nanmean(peBin(:,(indVal))), binMTdv.ccw(iSub,(indVal)), 1); 
    binMTdv.ccw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
    
    indVal = ~isnan(binMTdv.cw(iSub,:)); 
    pt = polyfit(nanmean(peBin(:,(indVal))), binMTdv.cw(iSub,(indVal)), 1); 
    binMTdv.cw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
end




set(figure(210),'position',[1 212 173 161]); clf; SP = subplot(1,1,1); hold on; 

hold on; 
plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv.ccw(:,ccwValInd)), nanstd(binMTdv.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv.ccw(:,ccwValInd))+nanstd(binMTdv.ccw(:,ccwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv.ccw(:,ccwValInd))-nanstd(binMTdv.ccw(:,ccwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv.ccw(:,ccwValInd)),'k-', 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(binMTdv.cw(:,cwValInd)), nanstd(binMTdv.cw(:,cwValInd),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv.cw(:,cwValInd))+nanstd(binMTdv.cw(:,cwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv.cw(:,cwValInd))-nanstd(binMTdv.cw(:,cwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv.cw(:,cwValInd)),'k-', 'color',cmap_c(3,:),'linewidth',1.5);
line([0 0],[-1 1],'linestyle','--','color','k'); 
line([-0.05 0.05],[0 0],'linestyle','--','color','k'); 
xlim([-0.05 0.05]); ylim([-.1 .1]); 
xlabel('Expectation'); ylabel('Decision variable'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-1 0 1], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')





asdfasdf
clear binMTdv_mimic
for iSub = 1:21
    for iBin = 1:(length(Rrange)-Bin_size)
        indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size)));
        if (sum(indCW)>=1)
            binMTdv_mimic.cw(iSub, iBin) = mean(MTdec_mimic{iSub}(indCW)); 
        else
            binMTdv_mimic.cw(iSub, iBin) = nan; 
        end
        if (sum(indCCW)>=1)
            binMTdv_mimic.ccw(iSub, iBin) = mean(MTdec_mimic{iSub}(indCCW)); 
        else
            binMTdv_mimic.ccw(iSub, iBin) = nan; 
        end
        binMTdv_mimic.all(iSub, iBin) = mean(MTdec_mimic{iSub}([indCW + indCCW]==1)); 
    end
end

cwValInd = sum(isnan(binMTdv_mimic.cw))<8; 
ccwValInd = sum(isnan(binMTdv_mimic.ccw))<8; 


% Substitute nan value
for iSub = 1:21
    indVal = ~isnan(binMTdv_mimic.ccw(iSub,:)); 
    pt = polyfit(nanmean(peBin(:,(indVal))), binMTdv_mimic.ccw(iSub,(indVal)), 1); 
    binMTdv_mimic.ccw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
    
    indVal = ~isnan(binMTdv_mimic.cw(iSub,:)); 
    pt = polyfit(nanmean(peBin(:,(indVal))), binMTdv_mimic.cw(iSub,(indVal)), 1); 
    binMTdv_mimic.cw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
end




set(figure(211),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 

hold on; 
plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv_mimic.ccw(:,ccwValInd)), nanstd(binMTdv_mimic.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv_mimic.ccw(:,ccwValInd))+nanstd(binMTdv_mimic.ccw(:,ccwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv_mimic.ccw(:,ccwValInd))-nanstd(binMTdv_mimic.ccw(:,ccwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv_mimic.ccw(:,ccwValInd)),'k-', 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(binMTdv_mimic.cw(:,cwValInd)), nanstd(binMTdv_mimic.cw(:,cwValInd),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv_mimic.cw(:,cwValInd))+nanstd(binMTdv_mimic.cw(:,cwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv_mimic.cw(:,cwValInd))-nanstd(binMTdv_mimic.cw(:,cwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv_mimic.cw(:,cwValInd)),'k-', 'color',cmap_c(3,:),'linewidth',1.5);
line([0 0],[-1 1],'linestyle','--','color','k'); 
% line([0 0],[-0.1 0.1],'linestyle','--','color','k'); 
line([-0.05 0.05],[0 0],'linestyle','--','color','k'); 
xlim([-0.05 0.05]); ylim([-1 1]); 
xlabel('Expectation'); ylabel('Decision variable'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-1 0 1], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')



%%

set(figure(212),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 

hold on; 
% plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv.ccw(:,ccwValInd)), nanstd(binMTdv.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_choiceStrong(1,:),0.5); 
% plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv.ccw(:,ccwValInd))+nanstd(binMTdv.ccw(:,ccwValInd),0,1)/sqrt(20),'k-', 'color',cmap_choiceStrong(1,:),'linewidth',0.7);
% plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv.ccw(:,ccwValInd))-nanstd(binMTdv.ccw(:,ccwValInd),0,1)/sqrt(20),'k-', 'color',cmap_choiceStrong(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv.ccw(:,ccwValInd)),'k-', 'color',cmap_choiceStrong(1,:),'linewidth',2);
% plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(binMTdv.cw(:,cwValInd)), nanstd(binMTdv.cw(:,cwValInd),0,1)/sqrt(20),cmap_choiceStrong(3,:),0.5); 
% plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv.cw(:,cwValInd))+nanstd(binMTdv.cw(:,cwValInd),0,1)/sqrt(20),'k-', 'color',cmap_choiceStrong(3,:),'linewidth',0.7);
% plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv.cw(:,cwValInd))-nanstd(binMTdv.cw(:,cwValInd),0,1)/sqrt(20),'k-', 'color',cmap_choiceStrong(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv.cw(:,cwValInd)),'k-', 'color',cmap_choiceStrong(3,:),'linewidth',2);


plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv_mimic.ccw(:,ccwValInd)), nanstd(binMTdv_mimic.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv_mimic.ccw(:,ccwValInd))+nanstd(binMTdv_mimic.ccw(:,ccwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv_mimic.ccw(:,ccwValInd))-nanstd(binMTdv_mimic.ccw(:,ccwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv_mimic.ccw(:,ccwValInd)),'k-', 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(binMTdv_mimic.cw(:,cwValInd)), nanstd(binMTdv_mimic.cw(:,cwValInd),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv_mimic.cw(:,cwValInd))+nanstd(binMTdv_mimic.cw(:,cwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv_mimic.cw(:,cwValInd))-nanstd(binMTdv_mimic.cw(:,cwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv_mimic.cw(:,cwValInd)),'k-', 'color',cmap_c(3,:),'linewidth',1.5);
line([0 0],[-1 1],'linestyle','--','color','k'); 
% line([0 0],[-0.1 0.1],'linestyle','--','color','k'); 
line([-0.05 0.05],[0 0],'linestyle','--','color','k'); 
xlim([-0.05 0.05]); ylim([-1 1]); 
xlabel('Expectation'); ylabel('Decision variable'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-1 0 1], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')


for iBin = 1:length(peBin) 
    ttestx(iBin,1) = ttest2(binMTdv_mimic.cw(:,iBin), binMTdv.cw(:,iBin)); 
    ttestx(iBin,2) = ttest2(binMTdv_mimic.ccw(:,iBin), binMTdv.ccw(:,iBin)); 
end

figure(); clf; hold on; 
plot(nanmean(binMTdv.cw), nanmean(binMTdv_mimic.cw), 'ko','markerfacecolor',cmap_c(3,:),'markeredgecolor','none','markersize',4); 
plot(nanmean(binMTdv.ccw), nanmean(binMTdv_mimic.ccw), 'ko','markerfacecolor',cmap_c(1,:),'markeredgecolor','none','markersize',4); 
line([-1 1],[-1 1],'linestyle','--','color','k');
% line([-1 1],[1 -1],'linestyle','--','color','k');
xlim([-1 1]); ylim([-1 1]);


[h,p,ci,stats] = ttest(nanmean(binMTdv_mimic.cw), nanmean(binMTdv.cw))
[h,p,ci,stats] = ttest(nanmean(binMTdv_mimic.ccw), nanmean(binMTdv.ccw))
[h,p,ci,stats] = ttest([nanmean(binMTdv_mimic.cw) -nanmean(binMTdv_mimic.ccw)], [nanmean(binMTdv.cw) -nanmean(binMTdv.ccw)])









%% OMPFC
Code this part again .... 







% Normalizing each bin's time-series
for iSub = 1:21
    DataMat_OMPFC{iSub}.sfm_norm = (DataMat_OMPFC{iSub}.sfm - repmat(mean(DataMat_OMPFC{iSub}.sfm,2), 1,500))./repmat(std(DataMat_OMPFC{iSub}.sfm,0,2), 1,500);
    DataMat_MT{iSub}.sfm_norm = (DataMat_MT{iSub}.sfm - repmat(mean(DataMat_MT{iSub}.sfm,2), 1,500))./repmat(std(DataMat_MT{iSub}.sfm,0,2), 1,500);
end

Bin_size = 20;
Bin_dt = 0.001;
peall = []; 
for iSub = 1:21
    peall = [peall PE_reg{iSub}];
end
maxPE = max(abs(peall)); 
Rrange = -maxPE:Bin_dt:maxPE; 

for iSub = 1:21
    for iT = 1:length(DataMat_OMPFC{iSub}.mimic_norm(1,:))
%         tp = polyfit(CorrVal.OMPFC_PE{iSub}(indK==3), DataMat_OMPFC{iSub}.sfm_norm((indK==3),iT), 1); 
        valV = ~isnan(CorrVal.OMPFC_PEmimic{iSub}); 
        tp = polyfit(CorrVal.OMPFC_PEmimic{iSub}(valV), DataMat_OMPFC{iSub}.mimic_norm((valV),iT), 1); 
        slopeDec_mimic{iSub}(1,iT) = tp(1); 
    end
    for iT = 1:length(DataMat_OMPFC{iSub}.sfm_norm(1,:))
%         tp = polyfit(CorrVal.OMPFC_PE{iSub}(indK==3), DataMat_OMPFC{iSub}.sfm_norm((indK==3),iT), 1); 
        valV = ~isnan(CorrVal.OMPFC_PE{iSub}); 
        tp = polyfit(CorrVal.OMPFC_PE{iSub}(valV), DataMat_OMPFC{iSub}.sfm_norm(valV,iT), 1); 
        slopeDec{iSub}(1,iT) = tp(1); 
    end
end

clear binBOLD peBin nBins
for iSub = 1:21
    for iBin = 1:(length(Rrange)-Bin_size)
        indCW = ((Dec_regorg{iSub} == 1) & (PE_regorg{iSub} > Rrange(iBin)) & (PE_regorg{iSub} <= Rrange(iBin+Bin_size))); 
        indCCW = ((Dec_regorg{iSub} ~= 1) & (PE_regorg{iSub} > Rrange(iBin)) & (PE_regorg{iSub} <= Rrange(iBin+Bin_size)));
        peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
        if (sum(indCW)>=1)
            binBOLD.cw(iSub, iBin) = mean(slopeDec{iSub}(indCW)); 
        else
            binBOLD.cw(iSub, iBin) = nan; 
        end
        if (sum(indCCW)>=1)
            binBOLD.ccw(iSub, iBin) = mean(slopeDec{iSub}(indCCW)); 
        else
            binBOLD.ccw(iSub, iBin) = nan; 
        end
        binBOLD.all(iSub, iBin) = mean(slopeDec{iSub}([indCW + indCCW]==1)); 
%         nBins(iSub,iBin,1) = sum(indCW); 
%         nBins(iSub,iBin,2) = sum(indCCW); 
    end
end

cwValInd = sum(isnan(binBOLD.cw))<2; 
ccwValInd = sum(isnan(binBOLD.ccw))<2; 

% % Substitute nan value
% for iSub = 1:21
%     indVal = ~isnan(binMTdv.ccw(iSub,:)); 
%     pt = polyfit(nanmean(peBin(:,(indVal))), binMTdv.ccw(iSub,(indVal)), 1); 
%     binMTdv.ccw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
%     
%     indVal = ~isnan(binMTdv.cw(iSub,:)); 
%     pt = polyfit(nanmean(peBin(:,(indVal))), binMTdv.cw(iSub,(indVal)), 1); 
%     binMTdv.cw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
% end

for iSub = 1:21
    for iBin = 1:(length(Rrange)-Bin_size)
        indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size)));
        peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
        if (sum(indCW)>=1)
            binBOLD_mimic.cw(iSub, iBin) = mean(slopeDec_mimic{iSub}(indCW)); 
        else
            binBOLD_mimic.cw(iSub, iBin) = nan; 
        end
        if (sum(indCCW)>=1)
            binBOLD_mimic.ccw(iSub, iBin) = mean(slopeDec_mimic{iSub}(indCCW)); 
        else
            binBOLD_mimic.ccw(iSub, iBin) = nan; 
        end
        binBOLD_mimic.all(iSub, iBin) = mean(slopeDec_mimic{iSub}([indCW + indCCW]==1)); 
%         nBins(iSub,iBin,1) = sum(indCW); 
%         nBins(iSub,iBin,2) = sum(indCCW); 
    end
end

cwValInd = sum(isnan(binBOLD_mimic.cw))<10; 
ccwValInd = sum(isnan(binBOLD_mimic.ccw))<10; 

% Substitute nan value
for iSub = 1:21
    indVal = ~isnan(binMTdv_mimic.ccw(iSub,:)); 
    pt = polyfit(nanmean(peBin(:,(indVal))), binMTdv_mimic.ccw(iSub,(indVal)), 1); 
    binMTdv_mimic.ccw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
    
    indVal = ~isnan(binMTdv_mimic.cw(iSub,:)); 
    pt = polyfit(nanmean(peBin(:,(indVal))), binMTdv_mimic.cw(iSub,(indVal)), 1); 
    binMTdv_mimic.cw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
end



set(figure(208),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 

pcw = polyfit(nanmean(peBin(:,cwValInd)), nanmean(binBOLD.cw(:,cwValInd)),1);
pccw = polyfit(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD.ccw(:,ccwValInd)),1);
pcw1 = polyfit(nanmean(peBin(:,cwValInd)), nanmean(binBOLD_mimic.cw(:,cwValInd)),1);
pccw1 = polyfit(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD_mimic.ccw(:,ccwValInd)),1);
% plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD_mimic.ccw(:,ccwValInd)), nanstd(binBOLD_mimic.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_c(1,:),0.5); 
% plot(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD_mimic.ccw(:,ccwValInd))+nanstd(binBOLD_mimic.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
% plot(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD_mimic.ccw(:,ccwValInd))-nanstd(binBOLD_mimic.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD.ccw(:,ccwValInd))-pccw(2), 'color',cmap_choiceStrong(1,:),'linewidth',2);
% plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(binBOLD_mimic.cw(:,cwValInd)), nanstd(binBOLD_mimic.cw(:,cwValInd),0,1)/sqrt(20),cmap_c(3,:),0.5); 
% plot(nanmean(peBin(:,cwValInd)), nanmean(binBOLD_mimic.cw(:,cwValInd))+nanstd(binBOLD_mimic.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
% plot(nanmean(peBin(:,cwValInd)), nanmean(binBOLD_mimic.cw(:,cwValInd))-nanstd(binBOLD_mimic.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binBOLD.cw(:,cwValInd))-pcw(2), 'color',cmap_choiceStrong(3,:),'linewidth',2);

plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD_mimic.ccw(:,ccwValInd))-pccw1(2), nanstd(binBOLD_mimic.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD_mimic.ccw(:,ccwValInd))-pccw1(2)+nanstd(binBOLD_mimic.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD_mimic.ccw(:,ccwValInd))-pccw1(2)-nanstd(binBOLD_mimic.ccw(:,ccwValInd),0,1)/sqrt(20), 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binBOLD_mimic.ccw(:,ccwValInd))-pccw1(2), 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(binBOLD_mimic.cw(:,cwValInd))-pcw1(2), nanstd(binBOLD_mimic.cw(:,cwValInd),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(binBOLD_mimic.cw(:,cwValInd))-pcw1(2)+nanstd(binBOLD_mimic.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binBOLD_mimic.cw(:,cwValInd))-pcw1(2)-nanstd(binBOLD_mimic.cw(:,cwValInd),0,1)/sqrt(20), 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binBOLD_mimic.cw(:,cwValInd))-pcw1(2), 'color',cmap_c(3,:),'linewidth',1.5);
line([0 0],[-1.5 1.5],'linestyle','--','color','k'); 
xlim([-0.05 0.05]); ylim([-1.5 1.5]); 
xlabel('Expectation'); ylabel('Weighted BOLD'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-1.5 0 1.5], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')
% 
% set(figure(209),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
% plot(nanmean(peBin(:,cwValInd)), mean(nBins(:,cwValInd,1),1),'color',cmap_c(3,:),'linewidth',1.3);
% plot(nanmean(peBin(:,ccwValInd)), mean(nBins(:,ccwValInd,2),1),'color',cmap_c(1,:),'linewidth',1.3);
% xlim([-0.05 0.05]); ylim([0 200]);
% set(SP, 'box', 'off', 'TickDir','out', 'YTick', [0 200], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')












































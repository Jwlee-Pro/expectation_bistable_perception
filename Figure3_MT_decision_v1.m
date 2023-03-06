clear all; 

addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2018_Mfin/Brain_residualAnalysis_2018June/PleaseWork/NeatVersion_library'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/Brain2Brain'); 

%% Load behavioral data 
Setting_param_bhv;
%% Generating the Model time-series (cHRF-convolved, raw)
Generating_model_TS; 

load('TS_OMPFC_all.mat'); 
load('TS_MT_all.mat'); 

PE_reg = cell(1,21); Dec_reg = cell(1,21); 
for iSub = 1:21
    for iS = 1:4
        PE_reg{iSub} = [PE_reg{iSub} mDataMat.PE_conv{iSub,iS}]; 
        Dec_reg{iSub} = [Dec_reg{iSub} mDataMat.Decision_conv{iSub,iS}];         
    end
    CorrVal.OMPFC_PE{iSub} = corr(DataMat_OMPFC{iSub}.sfm', PE_reg{iSub}');
    CorrVal.MT_D{iSub} = corr(DataMat_MT{iSub}.sfm', Dec_reg{iSub}');
end

%% Regressing out the U and H 
for iSub = 1:21
    for iS = 1:4
        for iVox = 1:length(DataMat_OMPFC{iSub}.sfm(:,1))
            TStemp = DataMat_OMPFC{iSub}.sfm(iVox,(1:125) + (iS-1)*125); 
            [b,bint,r,rint,stats] = regress(TStemp', [mDataMat.U_conv{iSub,iS}' ones(size(mDataMat.U_conv{iSub,iS}'))]); 
%             [b,bint,r,rint,stats] = regress(r, [mDataMat.H_conv{iSub,iS}' ones(size(mDataMat.H_conv{iSub,iS}'))]); 
            DataMat_OMPFC{iSub}.sfm_reg(iVox,(1:125) + (iS-1)*125) = r';
        end
    end
end

%% Normalizing each bin's time-series
for iSub = 1:21
    DataMat_OMPFC{iSub}.sfm_norm = (DataMat_OMPFC{iSub}.sfm_reg - repmat(mean(DataMat_OMPFC{iSub}.sfm_reg,2), 1,500))./repmat(std(DataMat_OMPFC{iSub}.sfm_reg,0,2), 1,500);
    DataMat_MT{iSub}.sfm_norm = (DataMat_MT{iSub}.sfm - repmat(mean(DataMat_MT{iSub}.sfm,2), 1,500))./repmat(std(DataMat_MT{iSub}.sfm,0,2), 1,500);
end



%% Sliding bin technique (pattern of activity difference between choice conditions) 
clear PatternSim PatternSim_rep
Bin_size = 15;
Bin_dt = 0.001;
peall = []; 
for iSub = 1:21
    peall = [peall PE_reg{iSub}];
end
maxPE = max(abs(peall)); 
Rrange = -maxPE:Bin_dt:maxPE; 
for iSub = 1:21
    for iBin = 1:(length(Rrange)-Bin_size)
        indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size)));
        
        peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
        if (sum(indCW)>1) & (sum(indCCW)>1)
            PatternSim.ompfc(iSub, iBin) = corr(mean(DataMat_OMPFC{iSub}.sfm_norm(:, indCW),2), mean(DataMat_OMPFC{iSub}.sfm_norm(:, indCCW),2)); 
            PatternSim.mt(iSub, iBin) = corr(mean(DataMat_MT{iSub}.sfm_norm(:, indCW),2), mean(DataMat_MT{iSub}.sfm_norm(:, indCCW),2)); 
        else
            PatternSim.ompfc(iSub, iBin) = nan; 
            PatternSim.mt(iSub, iBin) = nan; 
        end
    end
end


% Select the best voxels only
for iSub = 1:21
    [~, sI] = sort((CorrVal.OMPFC_PE{iSub})); 
    nOMPFC = length(CorrVal.OMPFC_PE{iSub}); 
    [~, sY] = sort((CorrVal.MT_D{iSub})); 
    nMT = length(CorrVal.MT_D{iSub}); 
    for iBin = 1:(length(Rrange)-Bin_size)
        indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size)));
        if (sum(indCW)>1) & (sum(indCCW)>1)
            PatternSim_rep.ompfc(iSub, iBin) = corr(mean(DataMat_OMPFC{iSub}.sfm_norm([sI(1:round(nOMPFC/10)); sI((nOMPFC-(round(nOMPFC/10))):nOMPFC)], indCW),2), mean(DataMat_OMPFC{iSub}.sfm_norm([sI(1:round(nOMPFC/10)); sI((nOMPFC-(round(nOMPFC/10))):nOMPFC)], indCCW),2)); 
            PatternSim_rep.mt(iSub, iBin) = corr(mean(DataMat_MT{iSub}.sfm_norm([sY(1:round(nMT/10)); sY((nMT-(round(nMT/10))):nMT)], indCW),2), mean(DataMat_MT{iSub}.sfm_norm([sY(1:round(nMT/10)); sY((nMT-(round(nMT/10))):nMT)], indCCW),2)); 
        else
            PatternSim_rep.ompfc(iSub, iBin) = nan; 
            PatternSim_rep.mt(iSub, iBin) = nan; 
        end
    end
end



figure(1); clf; 

SP = subplot(1,2,1); cla; hold on; 
plot(nanmean(peBin), nanmean(PatternSim.ompfc),'k.-'); 
plot(nanmean(peBin), nanmean(PatternSim_rep.ompfc),'r.-'); 
xlim([-0.025 0.025]); 
xlabel('Expectation'); ylabel('Pattern Similarity'); title('OMPFC'); 

SP = subplot(1,2,2); cla; hold on; 
plot(nanmean(peBin), nanmean(PatternSim.mt),'k.-'); 
plot(nanmean(peBin), nanmean(PatternSim_rep.mt),'r.-'); 
xlim([-0.025 0.025]); 
xlabel('Expectation'); ylabel('Pattern Similarity'); title('MT'); 




corrcrit = 0.134;
for iSub = 1:21
    corm(iSub,:) = sort(CorrVal.OMPFC_PE{iSub});
    sR(iSub) = sum(abs(CorrVal.OMPFC_PE{iSub}) > corrcrit);
end
[~,sortsR] = sort(sR); 

set(figure(204),'position',[4 699 282 106]); clf; SP = subplot(1,1,1); cla; hold on; 
imagesc(corm(sortsR,:)); colorbar; 
for rrSub = 1:21
    iSub = sortsR(rrSub); 
    [aa,ab] = sort(CorrVal.OMPFC_PE{iSub}); 
    line([0 0]+find(aa > corrcrit,1,'first'), [-0.5 0.5]+rrSub,'linestyle','-','color',[0 0 0]+0); 
%     line([find(aa > corrcrit,1,'first') 2537], [-0.5 -0.5]+rrSub,'linestyle','-','color',[0 0 0]+0.5); 
%     line([find(aa > corrcrit,1,'first') 2537], [0.5 0.5]+rrSub,'linestyle','-','color',[0 0 0]+0.5); 
    line([0 0]+find(aa < -corrcrit,1,'last'), [-0.5 0.5]+rrSub,'linestyle','-','color',[0 0 0]+0); 
%     line([0 find(aa < -corrcrit,1,'last')], [-0.5 -0.5]+rrSub,'linestyle','-','color',[0 0 0]+0.5); 
%     line([0 find(aa < -corrcrit,1,'last')], [0.5 0.5]+rrSub,'linestyle','-','color',[0 0 0]+0.5);
end
caxis([-0.5 0.5]);
xlim([0 length(CorrVal.OMPFC_PE{iSub})]); ylim([0.5 21.5]); 
load('/Users/joonwon/Desktop/Matlab_local/Brain_residualAnalysis_2018June/PleaseWork/cmap_choices'); 
colormap(mymap);
% xlabel('OMPFC voxels (sorted)'); %ylabel('Number of voxels'); 
set(SP, 'box', 'off', 'TickDir','out','XTick', [1 length(CorrVal.OMPFC_PE{iSub})],'FontSize', 10,'ycolor','none','xcolor','none')
set(SP, 'ydir','reverse'); 


% Weighted BOLD
clear decD
for iSub = 1:21
    Weights = CorrVal.OMPFC_PE{iSub};
    decD{iSub} = (Weights' * DataMat_OMPFC{iSub}.sfm_norm)/sqrt(sum(Weights.^2)); 
end
valNeural = [] ; 
valHD = [] ; 
for iS = 1:4
    valHD = [valHD [ones(1,123) 0 0]]; 
    valNeural = [valNeural [0 0 ones(1,123)]]; 
end
valNeural = valNeural==1;
valHD = valHD==1;
% Cross validation 
for iSub = 1:21
    PE_all_HDconsider{iSub}(valNeural) = PE_reg{iSub}(valHD); 
    for iS = 1:4
        PE_all_HDconsider{iSub}((1:2) + (iS-1)*125) = PE_reg{iSub}(3 + (iS-1)*125); 
    end
end
nT = 500; 
randInd = randperm(nT); 
nGG = 10; 
for iSub = 1:21
    clear Predc
    for iG = 1:nGG
        % Making model 
        indd = randInd((1:round(nT/nGG)) + (iG-1)*round(nT/nGG)); 
        orig_z = zeros(1,500); 
        orig_z(indd) = 1; 
        paramt = polyfit(PE_all_HDconsider{iSub}(orig_z==0), decD{iSub}(orig_z==0), 1); 
        Predc(orig_z==1) = (decD{iSub}(orig_z==1) - paramt(2))/paramt(1); 
    end
    Predicted{iSub} = Predc ;
end

clear cmv1 cmv2
for iSub = 1:21
    cmv1(iSub) = corr(PE_all_HDconsider{iSub}', decD{iSub}'); 
    cmv2(iSub) = corr(Predicted{iSub}', PE_all_HDconsider{iSub}'); 
end

set(figure(205),'position',[4 699 282 106]); clf; 
SP = subplot(1,1,1); cla; hold on; 
SP = subplot('position',[0.1300 0.1557 0.3 0.7693]); cla; hold on; 
bar(cmv2(sortsR),'facecolor',[0 0 0]+0.5,'edgecolor','none'); view(90,90);
errorbar(0.5, mean(cmv2), std(cmv2)/sqrt(20),'ko','markersize',4,'markerfacecolor','w'); 
plot(0.5, mean(cmv2),'ko', 'markerfacecolor','k','markeredgecolor','k','markersize',4)
xlim([0 22]); ylim([0 1]); 
set(SP, 'box', 'off', 'TickDir','out','XTick', [],'ytick',[],'FontSize', 10,'ycolor','k','xcolor','none')





%% Dissociate choice-explained variance and PE-explained variance
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
for iSub = 1:21
    for iVox = 1:length(DataMat_MT{iSub}.sfm(:,1))
        partCorrMat_mt.D{iSub}(iVox) = partcorr(DataMat_MT{iSub}.sfm(iVox,:)', Dec_reg{iSub}', PE_reg{iSub}'); 
        partCorrMat_mt.PE{iSub}(iVox) = partcorr(DataMat_MT{iSub}.sfm(iVox,:)', PE_reg{iSub}', Dec_reg{iSub}'); 
    end
    
    signPE = (CorrVal.MT_D{iSub}>0)*2 -1; 
    partCorrMat_mt.PE_signed{iSub} = partCorrMat_mt.PE{iSub}.*signPE'; 
    partCorrMat_mt.D_signed{iSub} = partCorrMat_mt.D{iSub}.*signPE'; 
end

diffVal = mean(partCorrMat.PE_signed) - mean(partCorrMat.D_signed); 
[~,sI] = sort(diffVal); 


diffM = [] ; 
for iSub = 1:21
    diffM = [diffM partCorrMat_mt.PE_signed{iSub}-partCorrMat_mt.D_signed{iSub}];
end


[aa1,a2]=ttest((partCorrMat.D_signed));
[aa2,a2]=ttest((partCorrMat.PE_signed));
for iVox = 1:length(aa1)
    if (aa1(iVox) == 1) & (aa2(iVox) == 1)
        indK(iVox) = 4; 
    elseif (aa1(iVox) == 1) & (aa2(iVox) ~= 1)
        indK(iVox) = 2; 
    elseif (aa1(iVox) ~= 1) & (aa2(iVox) == 1)
        indK(iVox) = 3; 
    else
        indK(iVox) = 1; % Both insignificant
    end
end


% cmap_type = [[0 0.5 0]+0.5;[0.5 0 0]+0.5; [0 0 0.5]+0.5;[0 0 0]]; 
cmap_type = [[0 0 0]+0.2;[0.5 0 0]+0.3; [0 0 0.5]+0.3;[0.5 0 0.5]+0.3]; 

set(figure(206),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
plot(mean(partCorrMat.PE_signed(:,indK==4)), mean(partCorrMat.D_signed(:,indK==4)),'kx','color',cmap_type(4,:),'markersize',3);
plot(mean(partCorrMat.PE_signed(:,indK==1)), mean(partCorrMat.D_signed(:,indK==1)),'kx','color',cmap_type(1,:),'markersize',3);
plot(mean(partCorrMat.PE_signed(:,indK==2)), mean(partCorrMat.D_signed(:,indK==2)),'kx','color',cmap_type(2,:),'markersize',3);
plot(mean(partCorrMat.PE_signed(:,indK==3)), mean(partCorrMat.D_signed(:,indK==3)),'kx','color',cmap_type(3,:),'markersize',3);
xlim([-0.05 0.1]); ylim([-0.05 0.1]); 
line([-0.05 0.1],[0 0],'linestyle','--','color','k'); 
line([0 0],[-0.05 0.1],'linestyle','--','color','k'); 
xlabel('Partcorr (PE)'); ylabel('Partcorr (D)'); 
[a1,a2]=corr(mean(partCorrMat.PE_signed)', mean(partCorrMat.D_signed)')
[a1,a2]=corr(mean(partCorrMat.PE)', mean(partCorrMat.D)')
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [0 .1], 'XTick', [0 .1], 'FontSize', 10,'color','none')

set(figure(210),'position',[1 447 80 161]); clf; SP = subplot(1,1,1); hold on; 
a = [sum(indK==1) sum(indK==2) sum(indK==3) sum(indK==4)]; 
h = bar([zeros(size(a)); a],'stacked'); 
for iCase = 1:4
    h(iCase).BarWidth = 0.4; 
    h(iCase).EdgeColor = 'none'; 
    h(iCase).FaceColor = cmap_type(iCase,:); 
end
xlim([1.5 2.5]); ylim([1 2537]); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [0 .1], 'XTick', [0 .1], 'FontSize', 10,'color','none','xcolor','none','ycolor','none')



%% Binned time-series


for iSub = 1:21
    for iT = 1:length(DataMat_OMPFC{iSub}.sfm_norm(1,:))
%         tp = polyfit(CorrVal.OMPFC_PE{iSub}(indK==3), DataMat_OMPFC{iSub}.sfm_norm((indK==3),iT), 1); 
        tp = polyfit(CorrVal.OMPFC_PE{iSub}, DataMat_OMPFC{iSub}.sfm_norm(:,iT), 1); 
        slopeDec{iSub}(1,iT) = tp(1); 
    end
end

Bin_size = 25;
Bin_dt = 0.001;
peall = []; 
for iSub = 1:21
    peall = [peall PE_reg{iSub}];
end
maxPE = max(abs(peall)); 
Rrange = -maxPE:Bin_dt:maxPE; 
for iSub = 1:21
    for iBin = 1:(length(Rrange)-Bin_size)
        indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size)));
        
        peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
        
        % Best voxel indices
        if (sum(indCW)>1) & (sum(indCCW)>1)
            PatternSim.ompfc(iSub, iBin) = corr(mean(DataMat_OMPFC{iSub}.sfm_norm(:, indCW),2), mean(DataMat_OMPFC{iSub}.sfm_norm(:, indCCW),2)); 
            PatternSim.mt(iSub, iBin) = corr(mean(DataMat_MT{iSub}.sfm_norm(:, indCW),2), mean(DataMat_MT{iSub}.sfm_norm(:, indCCW),2)); 
        else
            PatternSim.ompfc(iSub, iBin) = nan; 
            PatternSim.mt(iSub, iBin) = nan; 
        end
    end
end

clear binBOLD peBin nBins
numSelect = 100; 
for iSub = 1:21
%     PE_ind = find(indK == 3); 
%     [~, sI] = sort((CorrVal.OMPFC_PE{iSub}(PE_ind))); 
%     vPE_ind = PE_ind(sI); 
    for iBin = 1:(length(Rrange)-Bin_size)
        indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size)));
        peBin(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
        if (sum(indCW)>1)
            binBOLD.cw(iSub, iBin) = mean(slopeDec{iSub}(indCW)); 
        else
            binBOLD.cw(iSub, iBin) = nan; 
        end
        if (sum(indCCW)>1)
            binBOLD.ccw(iSub, iBin) = mean(slopeDec{iSub}(indCCW)); 
        else
            binBOLD.ccw(iSub, iBin) = nan; 
        end
        binBOLD.all(iSub, iBin) = mean(slopeDec{iSub}([indCW + indCCW]==1)); 
        nBins(iSub,iBin,1) = sum(indCW); 
        nBins(iSub,iBin,2) = sum(indCCW); 
    end
end


set(figure(207),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
pall = polyfit(nanmean(peBin), nanmean(binBOLD.all),1);
plot_patch_transparent(nanmean(peBin), nanmean(binBOLD.all), nanstd(binBOLD.all)/sqrt(20),cmap_c(2,:),0.5); 
plot(nanmean(peBin), nanmean(binBOLD.all)+nanstd(binBOLD.all,0,1)/sqrt(20), 'color',cmap_c(2,:),'linewidth',0.7);
plot(nanmean(peBin), nanmean(binBOLD.all)-nanstd(binBOLD.all,0,1)/sqrt(20), 'color',cmap_c(2,:),'linewidth',0.7);
plot(nanmean(peBin), nanmean(binBOLD.all), 'color',cmap_c(2,:),'linewidth',2);
xlim([-0.05 0.05]); ylim([-1.5 1.5]); 
xlabel('Expectation'); ylabel('Weighted BOLD'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-1.5 0 1.5], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')


set(figure(208),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
cwValInd = sum(isnan(binBOLD.cw))<2; 
ccwValInd = sum(isnan(binBOLD.ccw))<2; 
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
line([0 0],[-1.5 1.5],'linestyle','--','color','k'); 
xlim([-0.05 0.05]); ylim([-1.5 1.5]); 
xlabel('Expectation'); ylabel('Weighted BOLD'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-1.5 0 1.5], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')

set(figure(209),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); hold on; 
plot(nanmean(peBin(:,cwValInd)), mean(nBins(:,cwValInd,1),1),'color',cmap_c(3,:),'linewidth',1.3);
plot(nanmean(peBin(:,ccwValInd)), mean(nBins(:,ccwValInd,2),1),'color',cmap_c(1,:),'linewidth',1.3);
xlim([-0.05 0.05]); ylim([0 200]);
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [0 200], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')








% Load decoded decision variable from MT 
load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7/MTanalysis_deep/LeeLee_PriorEffect_workspace_PredCo20171127.mat','dec_final_fin'); 
MTdec = cell(1,21); 
for iSub = 1:21
    for iS = 1:4
        ta = dec_final_fin{4}{iSub,iS}; 
        MTdec{iSub} = [MTdec{iSub} log(ta)-log(1-ta)]; 
    end
end
clear binMTdv
for iSub = 1:21
    for iBin = 1:(length(Rrange)-Bin_size)
        indCW = ((Dec_reg{iSub} == 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size))); 
        indCCW = ((Dec_reg{iSub} ~= 1) & (PE_reg{iSub} > Rrange(iBin)) & (PE_reg{iSub} <= Rrange(iBin+Bin_size)));
        if (sum(indCW)>1)
            binMTdv.cw(iSub, iBin) = mean(MTdec{iSub}(indCW)); 
        else
            binMTdv.cw(iSub, iBin) = nan; 
        end
        if (sum(indCCW)>1)
            binMTdv.ccw(iSub, iBin) = mean(MTdec{iSub}(indCCW)); 
        else
            binMTdv.ccw(iSub, iBin) = nan; 
        end
        binMTdv.all(iSub, iBin) = mean(MTdec{iSub}([indCW + indCCW]==1)); 
    end
end



% Substitute nan value
for iSub = 1:21
    indVal = ~isnan(binMTdv.ccw(iSub,:)); 
    pt = polyfit(nanmean(peBin(:,(indVal))), binMTdv.ccw(iSub,(indVal)), 1); 
    binMTdv.ccw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
    
    indVal = ~isnan(binMTdv.cw(iSub,:)); 
    pt = polyfit(nanmean(peBin(:,(indVal))), binMTdv.cw(iSub,(indVal)), 1); 
    binMTdv.cw(iSub, (~indVal)) = pt(1)*peBin(iSub, (~indVal)) + pt(2);
end



% Loading the simulated data 
load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7/DataSimul_logPR'); 
RepSR = 15;  
Rrange = min(MatEvi_r{RepSR}(:,3)):Bin_dt:max(MatEvi_r{RepSR}(:,3));
clear mpeBin mBin_cw mBin_ccw mBin_rtcw mBin_rtccw
numCount = zeros((length(Rrange)-Bin_size),2); 
for iBin = 1:(length(Rrange)-Bin_size)
    mpeBin(iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
    indCW = ((MatEvi_r{RepSR}(:,1) == 1) & (MatEvi_r{RepSR}(:,3) > Rrange(iBin)) & (MatEvi_r{RepSR}(:,3) <= Rrange(iBin+Bin_size))); 
    indCCW = ((MatEvi_r{RepSR}(:,1) ~= 1) & (MatEvi_r{RepSR}(:,3) > Rrange(iBin)) & (MatEvi_r{RepSR}(:,3) <= Rrange(iBin+Bin_size)));
    
    TSa = MatEvi_r{RepSR}(:,2)'; 
    mBin_cw(iBin) = mean(TSa(indCW)); 
    mBin_ccw(iBin) = mean(TSa(indCCW)); 
    
    mBin_rtcw(iBin) = 1/mean(abs(TSa(indCW))); 
    mBin_rtccw(iBin) = 1/mean(abs(TSa(indCCW))); 
    
    numCount(iBin,1) = sum(indCW); 
    numCount(iBin,2) = sum(indCCW); 
end


set(figure(211),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); 
yyaxis right
% hold on; 
% PE = 0:0.01:1; 
% UU.cw = 1-PE; 
% UU.ccw = PE; 
% EU.cw = -PE.*log2(PE) - (1-PE).*log2(1-PE); 
% EU.ccw = -PE.*log2(PE) - (1-PE).*log2(1-PE); 
% % Model guide
% clear SPost
% xaxis = -10:0.001:11; 
% for iPE = 1:length(PE)
%     Posterior = normpdf(xaxis, PE(iPE)-0.5, 0.4)/sum(normpdf(xaxis, PE(iPE)-0.5, 0.4)); 
%     zeroInd = find(xaxis == 0); 
%     SPost.cw(iPE) = sum(xaxis(1:(zeroInd-1)).*Posterior(1:(zeroInd-1))); 
%     SPost.ccw(iPE) = sum(xaxis((zeroInd+1):length(Posterior)).*Posterior((zeroInd+1):length(Posterior))); 
% end
% plot(PE*0.12 - 0.06, SPost.cw,'k-','linestyle',':','color',cmap_choiceStrong(3,:),'linewidth',2);  
% plot(PE*0.12 - 0.06, SPost.ccw,'k-','linestyle',':','color',cmap_choiceStrong(1,:),'linewidth',2);
% ylim([-0.7 0.7]); 
% set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.1 0 0.1], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','none')

hold on; 
plot((mpeBin-0.5)/7, mBin_cw,'k-','linestyle',':','color',cmap_choiceStrong(1,:),'linewidth',2)
plot((mpeBin-0.5)/7, mBin_ccw,'k-','linestyle',':','color',cmap_choiceStrong(3,:),'linewidth',2)
ylim([-1.6 1.6])
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.1 0 0.1], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','none')

yyaxis left
hold on; 
plot_patch_transparent(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv.ccw(:,ccwValInd)), nanstd(binMTdv.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv.ccw(:,ccwValInd))+nanstd(binMTdv.ccw(:,ccwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv.ccw(:,ccwValInd))-nanstd(binMTdv.ccw(:,ccwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin(:,ccwValInd)), nanmean(binMTdv.ccw(:,ccwValInd)),'k-', 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin(:,cwValInd)), nanmean(binMTdv.cw(:,cwValInd)), nanstd(binMTdv.cw(:,cwValInd),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv.cw(:,cwValInd))+nanstd(binMTdv.cw(:,cwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv.cw(:,cwValInd))-nanstd(binMTdv.cw(:,cwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin(:,cwValInd)), nanmean(binMTdv.cw(:,cwValInd)),'k-', 'color',cmap_c(3,:),'linewidth',1.5);
% line([0 0],[0.35 0.65],'linestyle','--','color','k'); 
line([0 0],[-0.1 0.1],'linestyle','--','color','k'); 
line([-0.05 0.05],[0 0],'linestyle','--','color','k'); 
xlim([-0.05 0.05]); ylim([-0.1 0.1]); 
xlabel('Expectation'); ylabel('Decision variable'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.1 0 0.1], 'XTick', [-0.05 0 0.05], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')






%% Reaction time (RT) binning 
load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7/Data_allTS_15.mat'); 

Bin_size = 30;  
Bin_dt = 0.005; 
peall = []; 
for iSub = 1:21
    for iS = 1:4
        peall = [peall matSFM{iS,iSub}.result.matPrior];
    end
end
Rrange = min((peall)):Bin_dt:max((peall));

clear peBin_rt binRT 
nBin = 12; 
nDD = zeros((length(Rrange)-Bin_size),2); 
for iSub = 1:21
    load([session.ID{iSub} '_scanData_processed.mat'], 'paramSFM', 'paramMimic');
    pe_orig = []; d_orig = []; rt = [];
    for iS = 1:4
        paramSFM.track{iS} = Amb2twoAFC(paramSFM.track{iS}); 
        pe_orig = [pe_orig matSFM{iS,iSub}.result.matPrior];
        d_orig = [d_orig paramSFM.track{iS}(:,1)'];
        rt = [rt DataMat.reactionTime.zscore{iSub,iS}'];
    end
    
    for iBin = 1:(length(Rrange)-Bin_size)
        indCW = ((d_orig == 1) & (pe_orig > Rrange(iBin)) & (pe_orig <= Rrange(iBin+Bin_size))); 
        indCCW = ((d_orig ~= 1) & (pe_orig > Rrange(iBin)) & (pe_orig <= Rrange(iBin+Bin_size)));
        peBin_rt(iSub,iBin) = (Rrange(iBin)+Rrange(iBin+Bin_size))/2; 
        if (sum(indCW)>1)
            binRT.cw(iSub, iBin) = mean(rt(indCW)); 
        else
            binRT.cw(iSub, iBin) = nan; 
        end
        if (sum(indCCW)>1)
            binRT.ccw(iSub, iBin) = mean(rt(indCCW)); 
        else
            binRT.ccw(iSub, iBin) = nan; 
        end
        binRT.all(iSub, iBin) = mean(rt([indCW + indCCW]==1)); 
        nDD(iBin,1) = nDD(iBin,1) + sum(indCW); 
        nDD(iBin,2) = nDD(iBin,2) + sum(indCCW); 
    end
end

% Substitute nan value
for iSub = 1:21
    indVal = ~isnan(binRT.ccw(iSub,:)); 
    pt = polyfit(nanmean(peBin_rt(:,(indVal))), binRT.ccw(iSub,(indVal)), 1); 
    binRT.ccw(iSub, (~indVal)) = pt(1)*peBin_rt(iSub, (~indVal)) + pt(2);
    
    indVal = ~isnan(binRT.cw(iSub,:)); 
    pt = polyfit(nanmean(peBin_rt(:,(indVal))), binRT.cw(iSub,(indVal)), 1); 
    binRT.cw(iSub, (~indVal)) = pt(1)*peBin_rt(iSub, (~indVal)) + pt(2);
end
cwValInd = sum(isnan(binRT.cw))<10; 
ccwValInd = sum(isnan(binRT.ccw))<10; 


set(figure(212),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); 
% yyaxis right
% hold on; 
% plot(mpeBin, mBin_rtcw,'k-','linestyle',':','color',cmap_choiceStrong(1,:),'linewidth',2)
% plot(mpeBin, mBin_rtccw,'k-','linestyle',':','color',cmap_choiceStrong(3,:),'linewidth',2)
% xlim([0.2 0.8]);
% set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.3:0.3:0.6], 'XTick', [0.2:0.3:0.8], 'FontSize', 10,'color','none','xcolor','k','ycolor','none')


yyaxis right
hold on; 
plot(PE, 1./abs(SPost.cw),'k-','linestyle',':','color',cmap_choiceStrong(3,:),'linewidth',2);  
plot(PE, 1./abs(SPost.ccw),'k-','linestyle',':','color',cmap_choiceStrong(1,:),'linewidth',2);
% ylim([-0.7 0.7]); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.3:0.3:0.6], 'XTick', [0.2:0.3:0.8], 'FontSize', 10,'color','none','xcolor','k','ycolor','none')


yyaxis left
hold on; 
plot_patch_transparent(nanmean(peBin_rt(:,ccwValInd)), nanmean(binRT.ccw(:,ccwValInd)), nanstd(binRT.ccw(:,ccwValInd),0,1)/sqrt(20),cmap_c(1,:),0.5); 
plot(nanmean(peBin_rt(:,ccwValInd)), nanmean(binRT.ccw(:,ccwValInd))+nanstd(binRT.ccw(:,ccwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin_rt(:,ccwValInd)), nanmean(binRT.ccw(:,ccwValInd))-nanstd(binRT.ccw(:,ccwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(1,:),'linewidth',0.7);
plot(nanmean(peBin_rt(:,ccwValInd)), nanmean(binRT.ccw(:,ccwValInd)),'k-', 'color',cmap_c(1,:),'linewidth',1.5);
plot_patch_transparent(nanmean(peBin_rt(:,cwValInd)), nanmean(binRT.cw(:,cwValInd)), nanstd(binRT.cw(:,cwValInd),0,1)/sqrt(20),cmap_c(3,:),0.5); 
plot(nanmean(peBin_rt(:,cwValInd)), nanmean(binRT.cw(:,cwValInd))+nanstd(binRT.cw(:,cwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin_rt(:,cwValInd)), nanmean(binRT.cw(:,cwValInd))-nanstd(binRT.cw(:,cwValInd),0,1)/sqrt(20),'k-', 'color',cmap_c(3,:),'linewidth',0.7);
plot(nanmean(peBin_rt(:,cwValInd)), nanmean(binRT.cw(:,cwValInd)),'k-', 'color',cmap_c(3,:),'linewidth',1.5);
xlim([0.2 0.8]); ylim([-0.3 0.6]); 
xlabel('Expectation'); ylabel('Decision variable'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-0.3:0.3:0.6], 'XTick', [0.2:0.3:0.8], 'FontSize', 10,'color','none','xcolor','k','ycolor','k')



set(figure(213),'position',[1 447 173 161]); clf; SP = subplot(1,1,1); 
yyaxis left
hold on; 
vertex_x = [fliplr(nanmean(peBin_rt)) nanmean(peBin_rt)];
vertex_y = [zeros(size(nDD(:,1)'/21)) (nDD(:,1)'/21)  ];
v = [vertex_x; vertex_y]'; f = 1:length(v); 
patch('Faces',f,'Vertices',v,'facecolor',cmap_c(3,:),'edgecolor','none'); 
vertex_x = [fliplr(nanmean(peBin_rt)) nanmean(peBin_rt)];
vertex_y = [zeros(size(nDD(:,2)'/21)) -(nDD(:,2)'/21)  ];
v = [vertex_x; vertex_y]'; f = 1:length(v); 
patch('Faces',f,'Vertices',v,'facecolor',cmap_c(1,:),'edgecolor','none'); 
alpha(0.8)
xlabel('Expectation'); ylabel('Number of choices'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [-300 300], 'XTick', [0.2 0.5 0.8], 'FontSize', 10,'color','none','ycolor','k')
ylim([-300 300]); xlim([0.2 0.8]); 



yyaxis right 
cla; 
hold on; 
plot(mpeBin, numCount(:,1), 'k:','color', cmap_choiceStrong(1,:),'linewidth',2); 
plot(mpeBin, -numCount(:,2), 'k:','color', cmap_choiceStrong(3,:),'linewidth',2);
ylim([-300000 300000])
set(SP, 'box', 'off', 'TickDir','out', 'XTick', [0.2 0.5 0.8], 'FontSize', 10,'color','none','ycolor','none')



%% Predicting RT&choice from decoded DV
DurSet = [2.4 1.5 1.5 1.5 1.5 1.8 1.6 2.1 2.2 ...
        1.6 1.5 1.7 1.5 2.0 1.5 1.5 1.5 1.5 1.5 1.5 1.5]; 


nRange = 10 ; 
x_choice1 = nan(nRange, 21); 
y_choice1 = nan(nRange, 21); 
x_choice2 = nan(nRange, 21); 
y_choice2 = nan(nRange, 21); 


TR = 2.4; 
tmpTTR = (TR/2):TR:301; 
dv_TR = cell(21,4); 
HD = 6; 


MatResp_all = []; 
for iSub = 1:21
    load([session.ID{iSub} '_scanData_processed.mat'], 'paramSFM', 'paramMimic');
    pe_orig = []; d_orig = []; rt = []; MatResp = []; 
    for iS = 1:4
        paramSFM.track{iS} = Amb2twoAFC(paramSFM.track{iS}); 
        pe_orig = [pe_orig matSFM{iS,iSub}.result.matPrior];
        d_orig = [d_orig paramSFM.track{iS}(:,1)'];
        rt = [rt DataMat.reactionTime.zscore{iSub,iS}'];
        
        ta = dec_final_fin{4}{iSub,iS}; 
        dv_TR{iSub,iS} = log(ta)-log(1-ta); 
        
        % Resample DV values
        tmpT   = (DurSet(iSub)/2):DurSet(iSub):301; 
        
        resampledDV = interp1(tmpTTR, dv_TR{iSub,iS}, paramSFM.track{iS}(:,2)); 
        
        ttx = [resampledDV(round(HD/DurSet(iSub))+1:end) paramSFM.track{iS}(1:(end-round(HD/DurSet(iSub))),1) DataMat.reactionTime.zscore{iSub,iS}(1:(end-round(HD/DurSet(iSub))))]; 
%         ttx = [matSFM{iS,iSub}.result.matPrior(1:(end-round(HD/DurSet(iSub))))' paramSFM.track{iS}(1:(end-round(HD/DurSet(iSub))),1) DataMat.reactionTime.zscore{iSub,iS}(1:(end-round(HD/DurSet(iSub))))]; 
        MatResp = [MatResp; ttx]; 
    end
    Prediction(iSub) = (sum(MatResp(:,1)>0==MatResp(:,2)) + sum(-MatResp(:,1)<0==MatResp(:,2)))/length(MatResp(:,1)); 
    
    MatResp(MatResp(:,2)==-1,1) = - MatResp(MatResp(:,2)==-1,1) ;
    indval = ~isnan(sum(MatResp,2)); 
    cv(iSub) = corr(MatResp(indval,1), MatResp(indval,3)); 
    MatResp(MatResp(:,2)==-1,1) = - MatResp(MatResp(:,2)==-1,1) ;
    MatResp_all = [MatResp_all; MatResp]; 
    for ii = 1:nRange
        inds = find(MatResp(:,2)==1); 
        if ii == 1 
            id2 = MatResp(inds,1)<quantile(MatResp(inds,1), ii/nRange); 
        elseif ii == nRange
            id2 = MatResp(inds,1)>=quantile(MatResp(inds,1), (ii-1)/nRange); 
        else
            id2 = (MatResp(inds,1)>=quantile(MatResp(inds,1), (ii-1)/nRange)) & (MatResp(inds,1)<quantile(MatResp(inds,1), (ii)/nRange)); 
        end
        x_choice1(ii,iSub) = nanmean(MatResp(inds(id2),1)); 
        y_choice1(ii,iSub) = nanmean(MatResp(inds(id2),3)); 
        
        inds = find(MatResp(:,2)==-1); 
        if ii == 1 
            id2 = MatResp(inds,1)<quantile(MatResp(inds,1), ii/nRange); 
        elseif ii == nRange
            id2 = MatResp(inds,1)>=quantile(MatResp(inds,1), (ii-1)/nRange); 
        else
            id2 = (MatResp(inds,1)>=quantile(MatResp(inds,1), (ii-1)/nRange)) & (MatResp(inds,1)<quantile(MatResp(inds,1), (ii)/nRange)); 
        end
        x_choice2(ii,iSub) = nanmean(MatResp(inds(id2),1)) ; 
        y_choice2(ii,iSub) = nanmean(MatResp(inds(id2),3)) ; 
    end
end
% indval = ~isnan(sum(MatResp_all,2)) & MatResp_all(:,2)==1; 
% [coef, pval] = corr(MatResp_all(indval,1), MatResp_all(indval,3))
% 
% indval = ~isnan(sum(MatResp_all,2)) & MatResp_all(:,2)==-1; 
% [coef, pval] = corr(MatResp_all(indval,1), MatResp_all(indval,3))


figure(3); clf; 
SP = subplot(1,3,1);
hist(100*Prediction)
xlabel('DV-predicted Choice accuracy (%)'); 

SP = subplot(1,3,2);
hist(cv)
xlabel('Correlation betwee DV & RT'); 

[h,p,ci,stats] = ttest(cv)



SP = subplot(1,3,3); hold on; 
errorbar(nanmean(x_choice2'),nanmean(y_choice2'),std(y_choice2')/sqrt(20),'ko-','linewidth',2,'color',cmap_choiceStrong(3,:),'markerfacecolor','w','markersize',10); 
errorbar(nanmean(x_choice1'),nanmean(y_choice1'),std(y_choice1')/sqrt(20),'ko-','linewidth',2,'color',cmap_choiceStrong(1,:),'markerfacecolor','w','markersize',10); 
xlabel('Decoded DV');
ylabel('Reaction time (s)'); 



SP = subplot(1,3,3); cla; hold on; 
[~,sI] = sort(partc_DV); 
for isub = 1:(length(sub_list)-1)
%     if rawstatp{sI(isub)}.p<0.05
%         bar(isub, rawc_DV(sI(isub)),'Facecolor',[0 0 0]+0.3, 'Edgecolor','k'); 
%     else
%         bar(isub, rc_DV(sI(isub)),'Facecolor','w', 'Edgecolor','k'); 
%     end
    
    if statp{sI(isub)}.p<0.05
        bar(isub, partc_DV(sI(isub)),'Facecolor',[0 0 0]+0.3, 'Edgecolor','k'); 
    else
        bar(isub, partc_DV(sI(isub)),'Facecolor','w', 'Edgecolor','k'); 
    end
end
errorbar(52, nanmean(partc_DV), nanstd(partc_DV)/sqrt(length(sub_list)-1),'ko-','markerfacecolor','w'); 
view(90,90)
set(gca,'XAxisLocation','top','xdir','reverse'); 
xlabel('subject (sorted)');
ylabel('part correlation');
ylim([-0.35, 0.05])
xlim([-1, 55])






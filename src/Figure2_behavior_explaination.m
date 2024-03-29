%% Behavioral analysis & figure drawing

% Figure 2
% Correspondence of observed behavior and model

clear all; close all; clc; 


%% Parameters 
simParam.disparitySig = [0.4000 0.8000 1 1.3000 1.6000 2];
simParam.sampleSet = 1:50; 
simParam.nTrials = 200; 
simParam.nSession = 4000; 
simParam.alpha = 0.49:0.02:0.99; 


analysis_param.nSimPsession = 4; 
analysis_param.tFrag = 60; 
analysis_param.ShareFrag = linspace(0,analysis_param.tFrag,analysis_param.tFrag*2);
analysis_param.xT = 1:0.2:300; 

% Necessary library 
addpath('./library'); 
addpath('../data/behavior'); 

cmap_c=[102 204 255;204 204 204; 255 204 102]./255;
cmap_simul = [0.5 0 0]; 
cmap_repSub = [170 88 71]/255; 


% Loading general informations (cmap, session.ID)
ParamSetting_merge_fmri_pupil; 


idSub_rep = 8; 




%% Load simulated results (data not included in github repo.)
merge_do = 0; 

if merge_do == 1
    addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2018_Mf_Oct/Analysis_behavior/processed'); 
    addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2018_Mf_Oct/Analysis_behavior/processed0'); 
    for iDsp = 1:3
        for ialpha = 1:length(simParam.alpha)
            for iSR = 1:50
                if exist(['bistableDynamics_BistablePeak_alpha_' num2str(ialpha) '_SR_' num2str(iSR) '_DispSig_' num2str(iDsp) '.mat'], 'file') == 2
                    load(['bistableDynamics_BistablePeak_alpha_' num2str(ialpha) '_SR_' num2str(iSR) '_DispSig_' num2str(iDsp) '.mat']); 

                    simProcess_merge.MPR{iDsp,ialpha,iSR} = simProcess.MPR ;
                    simProcess_merge.RR{iDsp,ialpha,iSR} = simProcess.RR ;
                    simProcess_merge.avgDom(iDsp,ialpha,iSR) = simProcess.avgDom(2) ;
                    simProcess_merge.stdDom(iDsp,ialpha,iSR) = simProcess.stdDom(2) ;
                    simProcess_merge.RR_avg(iDsp,ialpha,iSR) = simProcess.RR_avg ;
                    simProcess_merge.SurvivalRate{iDsp,ialpha,iSR} = simProcess.SurvivalRate ;
                    simProcess_merge.CumTransitionRate{iDsp,ialpha,iSR} = simProcess.CumTransitionRate ;
                    simProcess_merge.SurvivalRate_all{1}(iDsp,ialpha,iSR,:) = simProcess.SurvivalRate_all{1} ;
                    simProcess_merge.SurvivalRate_all{2}(iDsp,ialpha,iSR,:) = simProcess.SurvivalRate_all{2} ;
                    simProcess_merge.SurvivalRate_all{3}(iDsp,ialpha,iSR,:) = simProcess.SurvivalRate_all{3} ;
                    simProcess_merge.CumTransitionRate_all{1}(iDsp,ialpha,iSR,:) = simProcess.CumTransitionRate_all{1} ;
                    simProcess_merge.CumTransitionRate_all{2}(iDsp,ialpha,iSR,:) = simProcess.CumTransitionRate_all{2} ;
                    simProcess_merge.CumTransitionRate_all{3}(iDsp,ialpha,iSR,:) = simProcess.CumTransitionRate_all{3} ;
                    simProcess_merge.CumProp(iDsp,ialpha,iSR,:) = simProcess.CumProp ;
                    simProcess_merge.GamFit(iDsp,ialpha,iSR,:) = simProcess.GamFit ;
                end
            end
        end
    end
    save('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2018_Mf_Oct/Analysis_behavior/bistableDynamics_BistablePeak_merged', 'simProcess_merge', '-v7.3'); 
else
    load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2018_Mf_Oct/Analysis_behavior/bistableDynamics_BistablePeak_merged'); 
end



%% Merging different sessions
clear beh_SFM beh_SFM_trialLength
ipupilSess = 12; 
beh_SFM = cell(1, length(session_merge.ID)); 
beh_SFM_trialLength = cell(1, length(session_merge.ID)); 
nSubject = length(session_merge.ID); 
for iSub = 1:nSubject
    
    if iSub < 13
        if length(session_merge.index{iSub}) < 3 
            for iSess = 1:length(session_merge.index{iSub})
                load([session_merge.index{iSub}{iSess} '_scanData_processed.mat'],'paramSFM');
                beh_SFM{iSub} = [beh_SFM{iSub} paramSFM.track];
                beh_SFM_trialLength{iSub}((1:length(paramSFM.track)) + length(beh_SFM{iSub})-4) = paramSFM.tDisplay + paramSFM.tBlank;
            end
        else
            load([session_merge.index{iSub} '_scanData_processed.mat'],'paramSFM');
            beh_SFM{iSub} = [beh_SFM{iSub} paramSFM.track];
            beh_SFM_trialLength{iSub}((1:length(paramSFM.track)) + length(beh_SFM{iSub})-4) = paramSFM.tDisplay + paramSFM.tBlank;
        end
        
    else
        for iSess = 1:ipupilSess
            try
                load([session_merge.index{iSub} '_scanData_' num2str(iSess) '_processed.mat'],'paramSFM');
                beh_SFM{iSub} = [beh_SFM{iSub} paramSFM.track];
                beh_SFM_trialLength{iSub}(iSess) = paramSFM.tDisplay + paramSFM.tBlank;
            catch
                fprintf(['missing sessions, iSub = ' num2str(iSub) ', Sess =' num2str(iSess) '\n']); 
            end
        end
    end
end


% Group summary
DataBehav.MPR = cell(1,nSubject);      % For dominance duration 
DataBehav.RR = cell(1,nSubject);       % Reversal rates 

clear paramSFM
for iSub = 1:nSubject
    
    for iS = 1:length(beh_SFM{iSub})
        beh_SFM{iSub}{iS} = Amb2twoAFC(beh_SFM{iSub}{iS});
        ValIndcw = (beh_SFM{iSub}{iS}(:,1)==1); 

        % Dominance duration & reversal rate 
        [a1, a2] = GenDominanceReversalRate(beh_SFM{iSub}{iS}(:,1));
        DataBehav.MPR{iSub} = [DataBehav.MPR{iSub}; a1*(beh_SFM_trialLength{iSub}(iS))]; 
        DataBehav.RR{iSub}  = [DataBehav.RR{iSub}; a2];
    end
    % Average dominance
    DataBehav.avgDom(iSub) = mean(abs(DataBehav.MPR{iSub}(:,2))); 
    DataBehav.stdDom(iSub) = std(abs(DataBehav.MPR{iSub}(:,2))); 
    
    % Average reversal rate 
    DataBehav.RR_avg(iSub) = mean(DataBehav.RR{iSub}(1:4)); 
    
    [xElapsed{iSub}, kb] = GenSurvivalRate(beh_SFM{iSub}, (nanmean(beh_SFM_trialLength{iSub})), analysis_param.tFrag);
    [~, kc] = GenCumTransitionRate_sess(beh_SFM{iSub}, (nanmean(beh_SFM_trialLength{iSub})), analysis_param.tFrag);
    
    % Survival & CumTransition raw 
    DataBehav.SurvivalRate{iSub}.merge(:,1) = kb{1};  % survival rate for cw percept
    DataBehav.SurvivalRate{iSub}.merge(:,2) = kb{2};  % survival rate for ccw percept
    DataBehav.SurvivalRate{iSub}.merge(:,3) = kb{3};  % survival rate for merged percept
    
    DataBehav.CumTransitionRate{iSub}.merge(:,1) = kc{1};  % cumulative transition rate for cw percept
    DataBehav.CumTransitionRate{iSub}.merge(:,2) = kc{2};  % cumulative transition rate for ccw percept
    DataBehav.CumTransitionRate{iSub}.merge(:,3) = kc{3};  % cumulative transition rate for merged percept
    
    
    % Survival & CumTransition resampled for comparison
    DataBehav.SurvivalRate_all{1}(iSub,:) = interp1(xElapsed{iSub}, kb{1}, analysis_param.ShareFrag); 
    DataBehav.SurvivalRate_all{2}(iSub,:) = interp1(xElapsed{iSub}, kb{2}, analysis_param.ShareFrag); 
    DataBehav.SurvivalRate_all{3}(iSub,:) = interp1(xElapsed{iSub}, kb{3}, analysis_param.ShareFrag); 
    
    DataBehav.CumTransitionRate_all{1}(iSub,:) = interp1(xElapsed{iSub}, kc{1}, analysis_param.ShareFrag); 
    DataBehav.CumTransitionRate_all{2}(iSub,:) = interp1(xElapsed{iSub}, kc{2}, analysis_param.ShareFrag); 
    DataBehav.CumTransitionRate_all{3}(iSub,:) = interp1(xElapsed{iSub}, kc{3}, analysis_param.ShareFrag); 
%     
    % Cumulative proportion (check for long tail) 
    for iTx = 1:length(analysis_param.xT)
        DataBehav.CumProp(iSub,iTx) = sum(abs(DataBehav.MPR{iSub}(:,2))<=analysis_param.xT(iTx))/length(DataBehav.MPR{iSub}(:,2)); 
    end
    [parmhat, parmci] = gamfit(abs(DataBehav.MPR{iSub}(:,2)));
    DataBehav.GamFit(iSub,:) = parmhat; 
end

%%
ExType = 2;
for iDsp = 1:3
    for ialpha = 1:length(simParam.alpha)
        for iSR = 1:length(simParam.sampleSet)

            % Individual fitting (vary all parameters)
            for iSub = 1:21
                % Calculate Hellinger distance 
                XX = gampdf(analysis_param.xT, DataBehav.GamFit(iSub,1), DataBehav.GamFit(iSub,2)); 
                XX_r = XX/sum(XX); 
                YY = gampdf(analysis_param.xT, simProcess_merge.GamFit(iDsp,ialpha,iSR,1), simProcess_merge.GamFit(iDsp,ialpha,iSR,2));
                YY_r = YY/sum(YY); 

                errorMerge_ind{iSub}(iDsp,ialpha, iSR) = sum((sqrt(XX_r) - sqrt(YY_r)).^2)/sqrt(2); % Hellinger distance
            end
        end
    end
end


% Find best parameter sets 
clear DataFit; 
for iType = ExType
    % Free = 3 
	for iSub = 1:21
        [~,I] = min(errorMerge_ind{iSub}(:));
        [DataFit.param_3dim(iSub,1), DataFit.param_3dim(iSub,2), DataFit.param_3dim(iSub,3)] = ind2sub(size(errorMerge_ind{iSub}),I);
    end
    % Free = 2 (Disp fixed to 2)
    for iSub = 1:21
        k1 = squeeze(errorMerge_ind{iSub}(2,:,:)); 
        [~,I] = min(k1(:));
        [DataFit.param_2dim(iSub,2), DataFit.param_2dim(iSub,3)] = ind2sub(size(k1),I);
        DataFit.param_2dim(iSub,1) = 2;
    end
    % Free = 1 (Disp fixed to 2; alpha fixed to 0.79)
    for iSub = 1:21
        k1 = squeeze(errorMerge_ind{iSub}(2,16,:)); 
        [~,I] = min(k1(:));
        [DataFit.param_1dim(iSub,3)] = ind2sub(size(k1),I);
        DataFit.param_1dim(iSub,1) = 2; DataFit.param_1dim(iSub,2) = 16;
    end
    % Free only alpha value 
    for iSub = 1:21
        k1 = squeeze(errorMerge_ind{iSub}(2,:,8:13)); 
        [~,I] = min(k1(:));
        [DataFit.param_alphaVary(iSub,2), DataFit.param_alphaVary(iSub,3)] = ind2sub(size(k1),I);
        DataFit.param_alphaVary(iSub,1) = 2; DataFit.param_alphaVary(iSub,3) = DataFit.param_alphaVary(iSub,3) + 7; 
    end
    for iSub = 1:21
        k1 = squeeze(errorMerge_ind{iSub}(2,:,12)); 
        [~,I] = min(k1(:));
        [DataFit.param_alphaVary1dim(iSub,2)] = ind2sub(size(k1),I);
        DataFit.param_alphaVary1dim(iSub,1) = 2; 
        DataFit.param_alphaVary1dim(iSub,3) = 12; 
    end
    
    % When fixating disp to other
    for iSub = 1:21
        k1 = squeeze(errorMerge_ind{iSub}(1,:,:)); 
        [~,I] = min(k1(:));
        [DataFit.param_2dim_1(iSub,2), DataFit.param_2dim_1(iSub,3)] = ind2sub(size(k1),I);
        DataFit.param_2dim_1(iSub,1) = 1;
        
        k1 = squeeze(errorMerge_ind{iSub}(2,:,:)); 
        [~,I] = min(k1(:));
        [DataFit.param_2dim_2(iSub,2), DataFit.param_2dim_2(iSub,3)] = ind2sub(size(k1),I);
        DataFit.param_2dim_2(iSub,1) = 2;
        
        k1 = squeeze(errorMerge_ind{iSub}(3,:,:)); 
        [~,I] = min(k1(:));
        [DataFit.param_2dim_3(iSub,2), DataFit.param_2dim_3(iSub,3)] = ind2sub(size(k1),I);
        DataFit.param_2dim_3(iSub,1) = 3;
    end
end

%%

Fig_size_cm{1} = [9 8]; 
Common_fontsize = 8;
Common_markersize = 9.5; 

BackgroundColor = [0 0 0]+1;
% BackgroundColor = [0 0 0]+1;
Fig_size_px{1} =round(cm2pixel(Fig_size_cm{1})); 

load('cmap_forFig2IndividualDiff_mild','cmap_ind'); 
[avgDom_sorted, IndSorted] = sort(DataBehav.avgDom); 




for iSub = 1:21
    BehPo.mean(iSub) = DataBehav.GamFit(iSub,1)*DataBehav.GamFit(iSub,2); 
    BehPo.std(iSub) = DataBehav.GamFit(iSub,1)*(DataBehav.GamFit(iSub,2)^2); 
    SimPo.mean_3dim(iSub) = simProcess_merge.GamFit(DataFit.param_3dim(iSub,1),DataFit.param_3dim(iSub,2),DataFit.param_3dim(iSub,3),1) * simProcess_merge.GamFit(DataFit.param_3dim(iSub,1),DataFit.param_3dim(iSub,2),DataFit.param_3dim(iSub,3),2); 
    SimPo.std_3dim(iSub) = simProcess_merge.GamFit(DataFit.param_3dim(iSub,1),DataFit.param_3dim(iSub,2),DataFit.param_3dim(iSub,3),1) * (simProcess_merge.GamFit(DataFit.param_3dim(iSub,1),DataFit.param_3dim(iSub,2),DataFit.param_3dim(iSub,3),2)^2); 
    SimPo.mean_2dim(iSub) = simProcess_merge.GamFit(DataFit.param_2dim(iSub,1),DataFit.param_2dim(iSub,2),DataFit.param_2dim(iSub,3),1) * simProcess_merge.GamFit(DataFit.param_2dim(iSub,1),DataFit.param_2dim(iSub,2),DataFit.param_2dim(iSub,3),2); 
    SimPo.std_2dim(iSub) = simProcess_merge.GamFit(DataFit.param_2dim(iSub,1),DataFit.param_2dim(iSub,2),DataFit.param_2dim(iSub,3),1) * (simProcess_merge.GamFit(DataFit.param_2dim(iSub,1),DataFit.param_2dim(iSub,2),DataFit.param_2dim(iSub,3),2)^2); 
    SimPo.mean_1dim(iSub) = simProcess_merge.GamFit(DataFit.param_1dim(iSub,1),DataFit.param_1dim(iSub,2),DataFit.param_1dim(iSub,3),1) * simProcess_merge.GamFit(DataFit.param_1dim(iSub,1),DataFit.param_1dim(iSub,2),DataFit.param_1dim(iSub,3),2); 
    SimPo.std_1dim(iSub) = simProcess_merge.GamFit(DataFit.param_1dim(iSub,1),DataFit.param_1dim(iSub,2),DataFit.param_1dim(iSub,3),1) * (simProcess_merge.GamFit(DataFit.param_1dim(iSub,1),DataFit.param_1dim(iSub,2),DataFit.param_1dim(iSub,3),2)^2); 
    
    SimPo.mean_alphaVary(iSub) = simProcess_merge.GamFit(DataFit.param_alphaVary(iSub,1),DataFit.param_alphaVary(iSub,2),DataFit.param_alphaVary(iSub,3),1) * simProcess_merge.GamFit(DataFit.param_alphaVary(iSub,1),DataFit.param_alphaVary(iSub,2),DataFit.param_alphaVary(iSub,3),2); 
    SimPo.std_alphaVary(iSub) = simProcess_merge.GamFit(DataFit.param_alphaVary(iSub,1),DataFit.param_alphaVary(iSub,2),DataFit.param_alphaVary(iSub,3),1) * (simProcess_merge.GamFit(DataFit.param_alphaVary(iSub,1),DataFit.param_alphaVary(iSub,2),DataFit.param_alphaVary(iSub,3),2)^2); 
    
end
BackgroundColor = 'w'; 


% Fit exponential functions
for iSub = 1:21
    % Survival probability 
    targetTS = DataBehav.SurvivalRate_all{3}(iSub,:)'; 
    IndVal = ~isnan(targetTS); 
    x = analysis_param.ShareFrag(IndVal)'; 
    g = fittype('0.5 + 0.5*exp(-a*x)');
    f0 = fit(x,targetTS(IndVal),g,'StartPoint',[1]);
    BehPo.SR(iSub) = f0.a; 

    targetTS = squeeze(simProcess_merge.SurvivalRate_all{3}(DataFit.param_1dim(iSub,1),DataFit.param_1dim(iSub,2),DataFit.param_1dim(iSub,3),:)); 
    IndVal = ~isnan(targetTS); 
    x = analysis_param.ShareFrag(IndVal)'; 
    g = fittype('0.5 + 0.5*exp(-a*x)');
    f0 = fit(x,targetTS(IndVal),g,'StartPoint',[1]);
    SimPo.SR_1dim(iSub) = f0.a; 
    
    targetTS = squeeze(simProcess_merge.SurvivalRate_all{3}(DataFit.param_alphaVary(iSub,1),DataFit.param_alphaVary(iSub,2),DataFit.param_alphaVary(iSub,3),:)); 
    IndVal = ~isnan(targetTS); 
    x = analysis_param.ShareFrag(IndVal)'; 
    g = fittype('0.5 + 0.5*exp(-a*x)');
    f0 = fit(x,targetTS(IndVal),g,'StartPoint',[1]);
    SimPo.SR_alphaVary(iSub) = f0.a; 
    
    targetTS = squeeze(simProcess_merge.SurvivalRate_all{3}(DataFit.param_3dim(iSub,1),DataFit.param_3dim(iSub,2),DataFit.param_3dim(iSub,3),:)); 
    IndVal = ~isnan(targetTS); 
    x = analysis_param.ShareFrag(IndVal)'; 
    g = fittype('0.5 + 0.5*exp(-a*x)');
    f0 = fit(x,targetTS(IndVal),g,'StartPoint',[1]);
    SimPo.SR_3dim(iSub) = f0.a; 
    
    targetTS = squeeze(simProcess_merge.SurvivalRate_all{3}(DataFit.param_2dim(iSub,1),DataFit.param_2dim(iSub,2),DataFit.param_2dim(iSub,3),:)); 
    IndVal = ~isnan(targetTS); 
    x = analysis_param.ShareFrag(IndVal)'; 
    g = fittype('0.5 + 0.5*exp(-a*x)');
    f0 = fit(x,targetTS(IndVal),g,'StartPoint',[1]);
    SimPo.SR_2dim(iSub) = f0.a; 

    % Cumulative transition probability 
    targetTS = DataBehav.CumTransitionRate_all{3}(iSub,:)'; 
    IndVal = ~isnan(targetTS); 
    x = analysis_param.ShareFrag(IndVal)'; 
    g = fittype('1-exp(-a*x)');
    f0 = fit(x,targetTS(IndVal),g,'StartPoint',[1]);
    BehPo.cTR(iSub) = f0.a; 


    targetTS = squeeze(simProcess_merge.CumTransitionRate_all{3}(DataFit.param_alphaVary(iSub,1),DataFit.param_alphaVary(iSub,2),DataFit.param_alphaVary(iSub,3),:)); 
    IndVal = ~isnan(targetTS); 
    x = analysis_param.ShareFrag(IndVal)'; 
    g = fittype('1-exp(-a*x)');
    f0 = fit(x,targetTS(IndVal),g,'StartPoint',[1]);
    SimPo.cTR_alphaVary(iSub) = f0.a; 
    
    targetTS = squeeze(simProcess_merge.CumTransitionRate_all{3}(DataFit.param_1dim(iSub,1),DataFit.param_1dim(iSub,2),DataFit.param_1dim(iSub,3),:)); 
    IndVal = ~isnan(targetTS); 
    x = analysis_param.ShareFrag(IndVal)'; 
    g = fittype('1-exp(-a*x)');
    f0 = fit(x,targetTS(IndVal),g,'StartPoint',[1]);
    SimPo.cTR_1dim(iSub) = f0.a; 
    
    targetTS = squeeze(simProcess_merge.CumTransitionRate_all{3}(DataFit.param_3dim(iSub,1),DataFit.param_3dim(iSub,2),DataFit.param_3dim(iSub,3),:)); 
    IndVal = ~isnan(targetTS); 
    x = analysis_param.ShareFrag(IndVal)'; 
    g = fittype('1-exp(-a*x)');
    f0 = fit(x,targetTS(IndVal),g,'StartPoint',[1]);
    SimPo.cTR_3dim(iSub) = f0.a; 
    
    targetTS = squeeze(simProcess_merge.CumTransitionRate_all{3}(DataFit.param_2dim(iSub,1),DataFit.param_2dim(iSub,2),DataFit.param_2dim(iSub,3),:)); 
    IndVal = ~isnan(targetTS); 
    x = analysis_param.ShareFrag(IndVal)'; 
    g = fittype('1-exp(-a*x)');
    f0 = fit(x,targetTS(IndVal),g,'StartPoint',[1]);
    SimPo.cTR_2dim(iSub) = f0.a; 
end





%% [Fig2.b] 
set(figure(1000),'position',[100 221 574 560]); 

SP=subplot('Position', [0.08     0.7205+0.05    0.13    0.14]); cla; hold on; 
for iSub = 1:nSubject
    bar(nSubject+1-iSub, avgDom_sorted(iSub), 'facecolor',cmap_ind(iSub,:),'edgecolor','none','linewidth',0.3); hold on; 
end
% line([0.5 21.5],[0 0]+1,'color','w','linewidth',0.75);
xlim([0 nSubject+1]); ylim([1 110]);  
set(gca,'YScale','log'); view(90,270); 
ylabel({'Mean duration (s)'}); xlabel('Subject index');
set(SP,'TickDir', 'out', 'box', 'off', 'FontSize', Common_fontsize, 'yTick',[1 5 10 30 100],'xTick',[],'xcolor','k')
set(gca,'color','w');


%% [Fig2.c] 
for iSub = 1:nSubject
    cmap_ind_sort(IndSorted(iSub),:) = cmap_ind(iSub,:); 
end
SP=subplot('Position', [0.29    0.7205+0.05    0.13    0.14]); cla; hold on;
plot(analysis_param.xT, DataBehav.CumProp,'color',[0 0 0]+0.7,'linewidth',1);
DataBehav.CumProp_sort = DataBehav.CumProp(IndSorted,:); 
for iSub = 1:length(session_merge.ID)
%     plot(analysis_param.xT, DataBehav.CumProp(iSub,:),'color',cmap_ind(iSub,:),'linewidth',1);
    plot(analysis_param.xT, DataBehav.CumProp(iSub,:),'color',cmap_ind_sort(iSub,:),'linewidth',1);
end
% plot(mean(DurSet(:,2))*8,1.03, 'marker','v','markeredgecolor','w','markerfacecolor','k','markersize',8); 
ylim([0 1.1]); xlim([1 max(analysis_param.xT)]); 
xlabel('Dom. duration (s)'); 
ylabel('Cumul. proportion'); 
set(gca,'XScale','log');
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [1 10 100 300], 'YTick', 0:1, 'FontSize', Common_fontsize);
set(gca,'color',BackgroundColor);
 
% plot(analysis_param.xT, DataBehav.CumProp(find(IndSorted==idSub_rep),:),'color',cmap_repSub,'linewidth',2);
% plot(analysis_param.xT, mean(DataBehav.CumProp),'color','k','linewidth',2.5);
% plot(analysis_param.xT, 0.5*ones(1,length(analysis_param.xT)),'linestyle','--','color','k'); 
% line([0 max(analysis_param.xT)], [0 0]+0.5,'linestyle','--','color','k'); 



% [Fig2.h] 
SP=subplot('Position', [0.35+0.21     0.7205+0.05    0.13    0.14]); cla; hold on;
for iSub = 1:length(session_merge.ID)
    plot(analysis_param.ShareFrag, DataBehav.SurvivalRate_all{3}(iSub,:),'color',cmap_ind_sort(iSub,:),'linewidth',1);
end
% plot(analysis_param.ShareFrag, DataBehav.SurvivalRate_all{3},'color',[0 0 0]+0.7,'linewidth',1);
% plot(analysis_param.ShareFrag, DataBehav.SurvivalRate_all{3}(find(IndSorted==idSub_rep),:),'color',cmap_repSub,'linewidth',2);
% plot(analysis_param.ShareFrag, mean(DataBehav.SurvivalRate_all{3}),'color','k','linewidth',2.5);
line([0 60], [0 0]+0.5,'linestyle','--','color','k'); 
ylim([0 1.1]); xlim([0 60]); 
xlabel('Elapsed time (s)'); 
ylabel('P_S'); 
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [0 20 40 60], 'YTick', 0:1, 'FontSize', Common_fontsize);
set(gca,'color',BackgroundColor);



% [Fig2.i] 
SP=subplot('Position', [0.35+0.42     0.7205+0.05    0.13    0.14]); cla; hold on;
for iSub = 1:length(session_merge.ID)
    plot(analysis_param.ShareFrag, DataBehav.CumTransitionRate_all{3}(iSub,:),'color',cmap_ind_sort(iSub,:),'linewidth',1);
end
% plot(analysis_param.ShareFrag, DataBehav.CumTransitionRate_all{3},'color',[0 0 0]+0.7,'linewidth',1);
% plot(analysis_param.ShareFrag, DataBehav.CumTransitionRate_all{3}(find(IndSorted==idSub_rep),:),'color',cmap_repSub,'linewidth',2);
% plot(analysis_param.ShareFrag, mean(DataBehav.CumTransitionRate_all{3}),'color','k','linewidth',2.5);
line([0 60], [0 0]+0.5,'linestyle','--','color','k'); 
ylim([0 1.1]); xlim([0 60]); 
xlabel('Elapsed time (s)'); 
ylabel('P_T'); 
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [0 20 40 60], 'YTick', 0:1, 'FontSize', Common_fontsize);
set(gca,'color',BackgroundColor);

drawnow; 




avgDomMat = squeeze(simProcess_merge.avgDom(2,:,5:18)); 

% [Fig2.c] 

SP=subplot('Position', [0.08    0.4+0.14    0.13    0.14]); cla; hold on;
% bone_rev = gray(80); %flipdim(gray(80),1); 
% imagesc(simParam.alpha, simParam.sampleSet(5:18),log(avgDomMat));
% caxis([0.9 4]); 
% colormap(bone_rev);
for iSub = 1:length(session_merge.ID)    
    plot(simParam.alpha(DataFit.param_alphaVary(IndSorted(iSub),2)), simParam.sampleSet(DataFit.param_alphaVary(IndSorted(iSub),3)),'marker','o','markerfacecolor',cmap_ind(iSub,:),'markeredgecolor','none', 'MarkerSize', 5)
end
xlim([0.4 1]); ylim([5 18]);  
ylabel('Sampling rate (\lambda)'); xlabel('Adapting rate (\alpha)'); 
set(SP,'TickDir', 'out', 'box', 'off', 'FontSize', Common_fontsize, 'XTick', [0.5 0.7 0.9 1], 'YTick', [5 10 15],'ycolor','k')






% [Fig2.e] 
Fix_lkhd = 2; 
clear DataSimul
for iSub = 1:nSubject
    XXX = simProcess_merge.MPR{Fix_lkhd, DataFit.param_alphaVary(iSub,2), DataFit.param_alphaVary(iSub,3)}/1.5*nanmean(beh_SFM_trialLength{iSub});
    for iTx = 1:length(analysis_param.xT)
        DataSimul.CumProp_sub{iSub}(iTx) = sum(abs(XXX(:,2))<=analysis_param.xT(iTx))/length(XXX); 
    end
end

SP=subplot('Position', [0.29    0.4+0.14    0.13    0.14]); cla; hold on;
for iSub = 1:length(session_merge.ID)
%     plot(analysis_param.xT, DataSimul.CumProp{gDataFit_dbm.param_ind(iSub,3)},'color',cmap_ind_sort(iSub,:),'linewidth',1);
    plot(analysis_param.xT, DataSimul.CumProp_sub{iSub},'color',cmap_ind_sort(iSub,:),'linewidth',1);
end
ylim([0 1.1]); xlim([1 max(analysis_param.xT)]); 
xlabel('Dom. duration (s)'); 
ylabel('Cumul. proportion'); 
set(gca,'XScale','log');
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [1 10 100 300], 'YTick', 0:1, 'FontSize', Common_fontsize);
set(gca,'color',BackgroundColor);



% [Fig2.j] 
SP=subplot('Position', [0.35+0.21     0.4+0.14   0.13    0.14]); cla; hold on;
for iSub = 1:length(session_merge.ID)
    plot(analysis_param.ShareFrag, squeeze(simProcess_merge.SurvivalRate_all{3}(Fix_lkhd,DataFit.param_alphaVary(iSub,2),DataFit.param_alphaVary(iSub,3),:)),'color',cmap_ind_sort(iSub,:),'linewidth',1);
end
line([0 60], [0 0]+0.5,'linestyle','--','color','k'); 
ylim([0 1.1]); xlim([0 60]); 
xlabel('Elapsed time (s)'); 
ylabel('P_S'); 
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [0 20 40 60], 'YTick', 0:1, 'FontSize', Common_fontsize);
set(gca,'color',BackgroundColor);

% [Fig2.k]
SP=subplot('Position', [0.35+0.42     0.4+0.14    0.13    0.14]); cla; hold on;
for iSub = 1:length(session_merge.ID) 
    plot(analysis_param.ShareFrag, squeeze(simProcess_merge.CumTransitionRate_all{3}(Fix_lkhd,DataFit.param_alphaVary(iSub,2),DataFit.param_alphaVary(iSub,3),:)),'color',cmap_ind_sort(iSub,:),'linewidth',1);
end
line([0 60], [0 0]+0.5,'linestyle','--','color','k'); 
ylim([0 1.1]); xlim([0 60]); 
xlabel('Elapsed time (s)'); 
ylabel('P_T'); 
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [0 20 40 60], 'YTick', 0:1, 'FontSize', Common_fontsize);
set(gca,'color',BackgroundColor);



% GOF
SP=subplot('Position', [0.08   0.2695    0.13    0.14]); cla; hold on;
for iSub = 1:nSubject 
    plot(BehPo.mean(iSub), SimPo.mean_alphaVary(iSub),'ko','markerfacecolor',cmap_ind_sort(iSub,:),'markeredgecolor','none','markersize',4);
end
line([3 60], [3 60],'linestyle','--','color','k'); 
set(gca,'XScale','log','YScale','log');
ylim([3 60]); xlim([3 60]); 
xlabel('Observed'); 
ylabel('Simulated'); 
title('Mean, \kappa\theta (s)','fontweight','bold');
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [3 10 60], 'YTick', [3 10 60], 'FontSize', Common_fontsize);
set(gca,'color',BackgroundColor);


SP=subplot('Position', [0.29    0.2695    0.13    0.14]); cla; hold on;
for iSub = 1:nSubject 
    plot(BehPo.std(iSub), SimPo.std_alphaVary(iSub),'ko','markerfacecolor',cmap_ind_sort(iSub,:),'markeredgecolor','none','markersize',4);
end
line([10 10000], [10 10000],'linestyle','--','color','k'); 
ylim([10 10000]); xlim([10 10000]); 
set(gca,'XScale','log','YScale','log');
xlabel('Observed'); 
ylabel('Simulated'); 
title('Variance, \kappa\theta^2 (s)','fontweight','bold');
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [10 10000], 'YTick', [10 10000], 'FontSize', Common_fontsize);
set(gca,'color',BackgroundColor); 



% [a1,a2] = corr(BehPo.std',SimPo.std')

SP=subplot('Position', [0.35+0.21    0.2695    0.13    0.14]); cla; hold on;
for iSub = 1:nSubject 
    plot(BehPo.SR(iSub), SimPo.SR_alphaVary(iSub),'ko','markerfacecolor',cmap_ind_sort(iSub,:),'markeredgecolor','none','markersize',4);
end
line([0.001 1], [0.001 1],'linestyle','--','color','k'); 
ylim([0.001 1]); xlim([0.001 1]); 
set(gca,'XScale','log','YScale','log');
xlabel('Observed'); 
ylabel('Simulated'); 
title('\tau_S (s)','fontweight','bold');
XTick = [0.001 0.01 0.1 1];
TickLabels{1} = '10^-^3'; TickLabels{2} = '10^-^2'; TickLabels{3} = '10^-^1'; TickLabels{4} = '1';
% set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [0.001 0.01 0.1 0.5], 'YTick', [0.001 0.01 0.1 0.5], 'FontSize', Common_fontsize);
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [0.001 0.01 0.1 1],'Xticklabel',TickLabels,'Yticklabel',TickLabels, 'YTick', [0.001 0.01 0.1 1], 'FontSize', Common_fontsize);
set(gca,'color',BackgroundColor); 


SP=subplot('Position', [0.35+0.42    0.2695    0.13    0.14]); cla; hold on;
for iSub = 1:nSubject 
    plot(BehPo.cTR(iSub), SimPo.cTR_alphaVary(iSub),'ko','markerfacecolor',cmap_ind_sort(iSub,:),'markeredgecolor','none','markersize',4);
end
line([0.001 1], [0.001 1],'linestyle','--','color','k'); 
ylim([0.001 1]); xlim([0.001 1]); 
set(gca,'XScale','log','YScale','log');
xlabel('Observed'); ylabel('Simulated'); 
title('\tau_T (s)','fontweight','bold');
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [0.001 0.01 0.1 1],'Xticklabel',TickLabels,'Yticklabel',TickLabels, 'YTick', [0.001 0.01 0.1 1], 'FontSize', Common_fontsize);
set(gca,'color',BackgroundColor); 

set(gcf, 'color','w');
FigHandle = gcf; 

save_figureToPDF(FigHandle, 'Figure2_BehaviorDynamics');


% SP=subplot('Position', [0.3500 0.7705 0.1300 0.1400]); hold on;
%     plot(analysis_param.xT, DataSimul.CumProp{20},'color','k','linewidth',2,'linestyle',':');
%     
% SP=subplot('Position', [0.5000 0.7705 0.1300 0.1400]); hold on;
%     plot(analysis_param.ShareFrag, squeeze(DataSimul.SurvivalRate_SIMmerge_rs(1,1,20,:)),'color','k','linewidth',2,'linestyle',':');
% 
% SP=subplot('Position', [0.7100 0.7705 0.1300 0.1400]); hold on;
%     plot(analysis_param.ShareFrag, squeeze(DataSimul.CumTransitionRate_SIMmerge_rs(1,1,20,:)),'color','k','linewidth',2,'linestyle',':');



%% Reaction time 
load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7/Data_allTS_15.mat'); 
load('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M7/MTanalysis_deep/LeeLee_PriorEffect_workspace_PredCo20171127.mat','matSFM','matMimic','matSFM_PredCo'); 

Bin_size = 15;
Bin_dt = 0.01;
peall = []; 
for iSub = 1:21
    for iS = 1:length(matSFM(:,iSub))
        peall = [peall matSFM{iS,iSub}.result.matPrior];
    end
end
Rrange = 0.3:Bin_dt:0.7; 


clear peBin_rt binRT coefval
nBin = 12; 
nDD = zeros(length(Rrange),2); 
pe_orig_all = []; 
d_orig_all = []; 
rt_all = []; 

for iSub = 1:21
    load([session.ID{iSub} '_scanData_processed.mat'], 'paramSFM', 'paramMimic');
    pe_orig = []; d_orig = []; rt = []; 
    for iS = 1:4
        paramSFM.track{iS} = Amb2twoAFC(paramSFM.track{iS}); 
        pe_orig = [pe_orig matSFM{iS,iSub}.result.matPrior];
        d_orig = [d_orig paramSFM.track{iS}(:,1)'];
        rt = [rt DataMat.reactionTime.zscore{iSub,iS}'];
    end
    pe_orig_all = [pe_orig_all pe_orig]; 
    d_orig_all = [d_orig_all d_orig]; 
    rt_all = [rt_all rt]; 
    
    [coefval(iSub), pval(iSub)] = corr((pe_orig'-0.5).*d_orig', rt','type','Spearman');
    
    
    for iBin = 1:length(Rrange)
        indCW = ((d_orig == 1) & (pe_orig > Rrange(iBin)-0.1) & (pe_orig <= Rrange(iBin)+0.1)); 
        indCCW = ((d_orig ~= 1) & (pe_orig > Rrange(iBin)-0.1) & (pe_orig <= Rrange(iBin)+0.1));
        peBin_rt(iSub,iBin) = Rrange(iBin);
        if (sum(indCW)>1)
            binRT.cw(iSub, iBin) = nanmean(rt(indCW)); 
        else
            binRT.cw(iSub, iBin) = nan; 
        end
        if (sum(indCCW)>1)
            binRT.ccw(iSub, iBin) = nanmean(rt(indCCW)); 
        else
            binRT.ccw(iSub, iBin) = nan; 
        end
        binRT.all(iSub, iBin) = mean(rt([indCW + indCCW]==1)); 
        nDD(iBin,1) = nDD(iBin,1) + sum(indCW); 
        nDD(iBin,2) = nDD(iBin,2) + sum(indCCW); 
    end
end

nanmean(coefval)
[h,p,ci,stats] = ttest(coefval)

[coefval_all, pval_all] = corr((pe_orig_all'-0.5).*d_orig_all', rt_all','type','Spearman')
[coefval_all, pval_all] = corr((pe_orig_all'-0.5).*d_orig_all', rt_all')
[coefval_all, pval_all] = corr((pe_orig_all'-0.5), d_orig_all','type','Spearman')
[coefval_all, pval_all] = corr((pe_orig_all'-0.5), d_orig_all')


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

figure(); hold on; 
plot(nanmean(peBin_rt), nanmean(binRT.ccw),'b-')
plot(nanmean(peBin_rt), nanmean(binRT.cw),'r-')

figure(40); clf; 
SP = subplot(1,2,1); cla; hold on; 
plot(pe_orig_all(d_orig_all==1), rt_all(d_orig_all==1),'k.'); 
plot(nanmean(peBin_rt), nanmean(binRT.cw),'w-','linewidth',3)
plot(nanmean(peBin_rt), nanmean(binRT.cw),'r-','linewidth',1.5)
SP = subplot(1,2,2); cla; hold on; 
plot(pe_orig_all(d_orig_all~=1), rt_all(d_orig_all~=1),'k.'); 
plot(nanmean(peBin_rt), nanmean(binRT.ccw),'w-','linewidth',3)
plot(nanmean(peBin_rt), nanmean(binRT.ccw),'b-','linewidth',1.5)


%% Figure 2ij (RT correspondence) 
figure(1001); clf; 
SP = subplot(1,1,1); cla; hold on; 
xval = [1-nanmean(peBin_rt) flip(1-nanmean(peBin_rt))]; 
temp = binRT.ccw; 
stdval  = nanstd(temp)/sqrt(length(session.ID)-1); 
meanval = nanmean(temp); 
yval = [-stdval+meanval +flip(stdval+meanval)] ; 
patch(xval, yval, cmap_c(1,:),'facealpha',0.4,'edgecolor','none'); 
plot(1-nanmean(peBin_rt),meanval,'color',cmap_c(1,:),'linewidth',1.5)

xval = [nanmean(peBin_rt) flip(nanmean(peBin_rt))]; 
temp = binRT.cw; 
stdval  = nanstd(temp)/sqrt(length(session.ID)-1); 
meanval = nanmean(temp); 
yval = [-stdval+meanval +flip(stdval+meanval)] ; 
patch(xval, yval, cmap_c(3,:),'facealpha',0.4,'edgecolor','none'); 
plot(nanmean(peBin_rt),meanval,'color',cmap_c(3,:),'linewidth',1.5)


set(gcf, 'color','w');

figure(); 

% [ss,sI] = sort(DataFit.param_alphaVary(:,2))
% [ss,sI] = sort(DataFit.param_alphaVary(:,3))
sI = IndSorted;

SP = subplot(1,2,1); cla; hold on; 
plot(nanmean(peBin_rt(sI(1:10),:)), nanmean(binRT.ccw(sI(1:10),:)),'b-')
plot(nanmean(peBin_rt(sI(11:end),:)), nanmean(binRT.ccw(sI(11:end),:)),'r-')
SP = subplot(1,2,2); cla; hold on; 
plot(nanmean(peBin_rt(sI(1:10),:)), nanmean(binRT.cw(sI(1:10),:)),'b-')
plot(nanmean(peBin_rt(sI(11:end),:)), nanmean(binRT.cw(sI(11:end),:)),'r-')
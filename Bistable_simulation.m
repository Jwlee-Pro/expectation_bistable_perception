%% Simulation movie 
%% Coded by Joonwon Lee 
%% Date 2017/6/15




clear all; 

% General Info
Common_fontsize = 12;
Common_markersize = 9; 
cmap_c=[102 204 255;204 204 204; 255 204 102]./255;
cmap_u=[207 207 207;255 111 207;128 250 0]./255;

% set(figure(2), 'Position', [100 100 994 875]); clf;
% set(gcf,'Color',[1 1 1]); 

addpath('/Volumes/Data_CSNL/people/JWL/Toolbox/NIfTI_toolbox/');
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Exp1_BOLD analysis/Behavior/'); % behavioral data
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/BayesianSampling/library/'); % parameter setting
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/PopulationCoding_MT/Core analysis codes'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Exp1_BOLD analysis/Library'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Code Library'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/PopulationCoding_MT'); 
addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/SfM_FinalAnalysis_version/Tuned_fin_Nov'); 


%% Model parameters 

    % [1] Likelihoods 
    LKHDparam.nVox = 1; 
    LKHDparam.axisVal = -10:0.1:10; 
%     LKHDparam.Voxmean = rand(2,LKHDparam.nVox)*2 - 1; % 1st row: cw, 2nd: ccw
    LKHDparam.Voxmean = [1; -1]; % 1st row: cw, 2nd: ccw
    LKHDparam.Voxsigma = (rand(1,LKHDparam.nVox)+1)/2; % Assume same std for cw- & ccw-preferring voxels 
    LKHDparam.Voxsigma = 1.5; % Assume same std for cw- & ccw-preferring voxels 
    for iVox = 1:LKHDparam.nVox
        LKHDparam.pdf_cw(iVox,:)  = normpdf(LKHDparam.axisVal, LKHDparam.Voxmean(1,iVox), LKHDparam.Voxsigma(iVox))./sum(normpdf(LKHDparam.axisVal, LKHDparam.Voxmean(1,iVox), LKHDparam.Voxsigma(iVox)));
        LKHDparam.pdf_ccw(iVox,:) = normpdf(LKHDparam.axisVal, LKHDparam.Voxmean(2,iVox), LKHDparam.Voxsigma(iVox))./sum(normpdf(LKHDparam.axisVal, LKHDparam.Voxmean(2,iVox), LKHDparam.Voxsigma(iVox)));
        LKHDparam.pdf_all(iVox,:) = (LKHDparam.pdf_cw(iVox,:) + LKHDparam.pdf_ccw(iVox,:))/2;
    end

    % [2] Prior 
    PRparam.alpha = 0.8; 
    PRparam.vectPi = 0:0.01:1 ;
    PRparam.vectPrior_init = ones(size(PRparam.vectPi))./length(PRparam.vectPi);

    % [3] Posterior & sampled posterior
    POparam.sampleSet = 1:50; 
    POparam.sampleRep = 20; 
    POparam.DecisionBound = 7; 
    

    % [4] Simulation parameters 
    SimParam.nTrials = 10000; 

    
%% Model simulation "main"
iVox = 1; 
clear RespMat result; 

for iN = POparam.sampleRep
    for iT = 1:SimParam.nTrials
        % (1) Prior calculation 
        if iT == 1
            matPrior(iT,:) = PRparam.vectPrior_init;
        else
            matPrior(iT,:)=((matPOST(iT-1,:)*PRparam.alpha) + (1-PRparam.alpha)*PRparam.vectPrior_init(locMax(iT-1)))/sum((matPOST(iT-1,:)*PRparam.alpha) + (1-PRparam.alpha)*PRparam.vectPrior_init(locMax(iT-1)));
        end
        result.matPrior(iT)=sum(PRparam.vectPi.*matPrior(iT,:));
        
        % (2) Random value extraction 
        vectRand = rand(POparam.sampleSet(iN),1);
        tL_cw    = find(vectRand<result.matPrior(iT));  
        tL_ccw   = find(vectRand>=result.matPrior(iT));
        inDpr    = vectRand<result.matPrior(iT); 

        vectOutcome_A=[];
        if(~isempty(tL_cw))
            vectOutcome_A(tL_cw) = randn(size(tL_cw)).*LKHDparam.Voxsigma + LKHDparam.Voxmean(1);
        end
        if(~isempty(tL_ccw))
            vectOutcome_A(tL_ccw) = randn(size(tL_ccw)).*LKHDparam.Voxsigma + LKHDparam.Voxmean(2);
        end
        RespMat{iN}(iT,1) = mean(vectOutcome_A);
        RespMat_allsample{iN}(iT,:) = (vectOutcome_A);
        RespMat_whichDist{iN}(iT,:) = (inDpr);
        
        SumLogPR = 0; SumMat = [0];
        
        for iSample = 1:POparam.sampleSet(iN)
          
            % LogPR (within-trial) 
            tempLoc = find(LKHDparam.axisVal <= vectOutcome_A(iSample),1,'last');
            valIndset = [tempLoc tempLoc+1] ;
            % Posterior calculation 
            if valIndset(1) == length(LKHDparam.pdf_cw)
                val1_cw = LKHDparam.pdf_cw(iVox,valIndset(1)); 
                val1_ccw = LKHDparam.pdf_ccw(iVox,valIndset(1)); 
            else
                val1_cw = (LKHDparam.pdf_cw(iVox,valIndset(1)) + LKHDparam.pdf_cw(iVox,valIndset(2)))/2; 
                val1_ccw = (LKHDparam.pdf_ccw(iVox,valIndset(1)) + LKHDparam.pdf_ccw(iVox,valIndset(2)))/2; 
            end
            
            SumLogPR = SumLogPR + log(val1_cw/val1_ccw); 
            SumMat = [SumMat SumLogPR]; 
            % Save data 
            RespMat_logPR{iN}(iT,iSample) = log(val1_cw/val1_ccw); 
            RespMat_logPR_cum{iN}(iT,iSample) = SumLogPR; 
        end


        % Saving the data
        RespMat{iN}(iT,2) = SumMat(end)>0; 
%         RespMat_logPR{iN}(iT,:) = (inDpr);
        
        % (3) Update prior depending on the decision
        matLKHD(iT,:)=(PRparam.vectPi.^sum(RespMat{iN}(iT,2))).*((1-PRparam.vectPi).^(1-sum(RespMat{iN}(iT,2))));
        matPOST(iT,:)=(matLKHD(iT,:).*matPrior(iT,:))./sum(matLKHD(iT,:).*matPrior(iT,:));
        result.matPOST(iT,:)=matPOST(iT,:);
        
        locMax(iT)=find(matPOST(iT,:)==max(matPOST(iT,:)),1,'first');
        result.maxPi(iT) = sum(PRparam.vectPi.*matPrior(iT,:)); % DBM (general method), not MAP
        result.maxPi_map(iT) = PRparam.vectPi(locMax(iT)); % MAP method 
                
        
        % (4) Obtain the uncertainty signals 
        result.UU(iT) =abs(result.maxPi(iT) - RespMat{iN}(iT,2)); 
        result.EU(iT) = -result.maxPi(iT)*log2(result.maxPi(iT)) - ((1-result.maxPi(iT))*log2(1-result.maxPi(iT)));

    end
end
                
% save('SimulationForResult_50000')

RespMat{iN}((RespMat{iN}(:,2) == 0),2) = -1;  

% extract run length statistics
matRun=[];
cntRun=0;
% 1st trila
currentCh=RespMat{iN}(1,2);
cntRun=cntRun+1;
cntTrial=1;
flagChange=0;
lengthRun=1;
while(cntTrial<length(RespMat{iN}(:,1)))
    while(flagChange==0)
        cntTrial=cntTrial+1;
        if(RespMat{iN}(cntTrial,2)==currentCh)

           lengthRun=lengthRun+1;
           if cntTrial == length(RespMat{iN}(:,1))
               matRun(cntRun,1)=currentCh;
               matRun(cntRun,2)=lengthRun;
                break; 
           end
        else

            matRun(cntRun,1)=currentCh;
            matRun(cntRun,2)=lengthRun;

            cntRun=cntRun+1;
            lengthRun;
            lengthRun=1;
            currentCh=RespMat{iN}(cntTrial,2);
            break;
        end
        if(cntTrial==length(RespMat{iN}(:,1)))
            break;
        end
    end

    flagChange=0;
end


lengthStable = 8 ; 
% sliding window
vectStable = zeros(1,length(RespMat{iN}(:,1))); 
for iT = 1:(length(RespMat{iN}(:,1))-lengthStable+1)
    indRange = iT:iT+lengthStable-1; 
    if (sum(RespMat{iN}(indRange,2)) == lengthStable) || (sum(RespMat{iN}(indRange,2)) == -lengthStable)
        vectStable(indRange) = 1; 
    end
end
% Detect single change
for iT = 1:(length(vectStable)-1)
    if (vectStable(iT) == 1) && (RespMat{iN}(iT,2)~=RespMat{iN}(iT+1,2))
        vectStable(iT+1) = 0 ; 
    end
end

        

matRun(:,3)=cumsum(matRun(:,2)); %cumulative trial end
matRun(1,4)=1;
matRun(2:end,4)=matRun(1:end-1,3)+1; % cumuylative trial start




% Calculate RT 
for iT = 1:length(RespMat_logPR_cum{iN}(:,1))
    if isempty(find(abs(RespMat_logPR_cum{iN}(iT,:)) > POparam.DecisionBound,1,'first'))
       RT_model(iT) = length(RespMat_logPR_cum{iN}(iT,:)); 
       RT_model_plot(iT) = length(RespMat_logPR_cum{iN}(iT,:)) + rand(1)*2 - 0.5; 
    else
       RT_model(iT) = find(abs(RespMat_logPR_cum{iN}(iT,:)) > POparam.DecisionBound,1,'first'); 
       RT_model_plot(iT) = find(abs(RespMat_logPR_cum{iN}(iT,:)) > POparam.DecisionBound,1,'first'); 
    end
end


asdf

%% Draw figure for simple insets
set(figure(3), 'Position', [2 496 541 298]); clf;
SP=subplot(2,3,1); cla; hold on; 
LKHD_sense = 0.5*LKHDparam.pdf_cw + 0.5*LKHDparam.pdf_ccw; 
plot(LKHDparam.axisVal, LKHD_sense, 'color',cmap_c(2,:),'linewidth',3); hold on; 
% plot(LKHDparam.axisVal, 0.8*0.1*LKHDparam.pdf_ccw,'color',cmap_c(3,:),'linewidth',3);
% line([0 0],[0 0.034],'linestyle','--','color','k'); 
xlim([-7 7]); ylim([-0.0007 0.03]); 
title('p(S|F)'); 
set(SP, 'box', 'off', 'TickDir', 'out','XTick',[0],'YTick',[],'FontSize', Common_fontsize,'ycolor','w');


SP=subplot(2,3,2); cla; hold on; 
plot(LKHDparam.axisVal, (LKHDparam.axisVal-min(LKHDparam.axisVal)).^2, 'color',cmap_c(3,:),'linewidth',3); hold on; 
plot(LKHDparam.axisVal, (-LKHDparam.axisVal-min(-LKHDparam.axisVal)).^2, 'color',cmap_c(1,:),'linewidth',3); hold on; 
% plot(LKHDparam.axisVal, 0.8*0.1*LKHDparam.pdf_ccw,'color',cmap_c(3,:),'linewidth',3);
xlim([-10 10]); ylim([-10 420]); 
set(SP, 'box', 'off', 'TickDir', 'out','XTick',[0],'YTick',[],'FontSize', Common_fontsize,'ycolor','w');
title('p(F|\Omega)'); 

SP=subplot(2,3,3); cla; hold on; 
plot(LKHDparam.axisVal, LKHD_sense.*((LKHDparam.axisVal-min(LKHDparam.axisVal)).^2), 'color',cmap_c(3,:),'linewidth',3); hold on; 
plot(LKHDparam.axisVal, ((-LKHDparam.axisVal-min(-LKHDparam.axisVal)).^2).*LKHD_sense, 'color',cmap_c(1,:),'linewidth',3); hold on; 
xlim([-7 7]); ylim([-0.1 3]); 
set(SP, 'box', 'off', 'TickDir', 'out','XTick',[0],'YTick',[],'FontSize', Common_fontsize,'ycolor','w');
title('p(S|\Omega)'); 


SP=subplot(2,3,5); cla; hold on; 
plot(LKHDparam.axisVal, ((LKHDparam.axisVal-min(LKHDparam.axisVal)).^2)*0.8, 'color',cmap_c(3,:),'linewidth',3); hold on; 
plot(LKHDparam.axisVal, ((-LKHDparam.axisVal-min(-LKHDparam.axisVal)).^2)*0.2, 'color',cmap_c(1,:),'linewidth',3); hold on; 
plot(LKHDparam.axisVal, ((LKHDparam.axisVal-min(LKHDparam.axisVal)).^2)*0.8 + ((-LKHDparam.axisVal-min(-LKHDparam.axisVal)).^2)*0.2, 'color','k','linewidth',3); hold on; 
xlim([-10 10]); ylim([-10 420]); 
set(SP, 'box', 'off', 'TickDir', 'out','XTick',[0],'YTick',[],'FontSize', Common_fontsize,'ycolor','w');
title('p(\Omega|F) when p(\Omega=cw)=0.8'); 


SP=subplot(2,3,6); cla; hold on; 
plot(LKHDparam.axisVal, LKHD_sense.*((LKHDparam.axisVal-min(LKHDparam.axisVal)).^2)*0.8 + LKHD_sense.*((-LKHDparam.axisVal-min(-LKHDparam.axisVal)).^2)*0.2, 'color','k','linewidth',3); hold on; 
RR = [max(LKHD_sense.*((LKHDparam.axisVal-min(LKHDparam.axisVal)).^2)*0.8 + LKHD_sense.*((-LKHDparam.axisVal-min(-LKHDparam.axisVal)).^2)*0.2) min(LKHD_sense.*((LKHDparam.axisVal-min(LKHDparam.axisVal)).^2)*0.8 + LKHD_sense.*((-LKHDparam.axisVal-min(-LKHDparam.axisVal)).^2)*0.2)];
xlim([-7 7]); %ylim([0 420]); 
% line([0 0],[0 3],'linestyle','--','color','k'); 
set(SP, 'box', 'off', 'TickDir', 'out','XTick',[0],'YTick',[],'FontSize', Common_fontsize,'ycolor','w');
title('p(\Omega|S)'); 

% Random value extraction (sampled posterior)
vectRand = rand(30,1);
tL_cw    = find(vectRand<0.8);  
tL_ccw   = find(vectRand>=0.8);
vectOutcome_A=[];
if(~isempty(tL_cw))
    vectOutcome_A(tL_cw) = randn(size(tL_cw)).*LKHDparam.Voxsigma + LKHDparam.Voxmean(1);
end
if(~isempty(tL_ccw))
    vectOutcome_A(tL_ccw) = randn(size(tL_ccw)).*LKHDparam.Voxsigma + LKHDparam.Voxmean(2);
end
hold on; 
h1 = hist(vectOutcome_A,-7:1:7);
h2 = bar(-7:1:7, RR(1)*h1/max(h1),'Facecolor',[0 0 0]+0.5, 'Edgecolor','none','Facealpha',0.6);
ylim([0 3]); 

%% Draw Figure 
figure(102); clf; 
FigureNumber = 102; 
Fig_size_cm{1} = [11 14]; 
Common_fontsize = 12;
Common_markersize = 9.5; 
Fig_size_px{1} =round(cm2pixel(Fig_size_cm{1})); 
set(figure(FigureNumber), 'Position', [1 409 3*Fig_size_px{1}]); clf;
set(gcf,'Color',[1 1 1]); 

SP = subplot('position', [0.150    0.88    0.70    0.09]); cla; hold on;
image(1:100,linspace(1,-1,101),flipud(255*(1-result.matPOST'./max(max(result.matPOST))))); hold on;
colormap(gray(256));
line([0 101],[0 0],'linestyle','--','color','k'); 
plot(1:100,result.maxPi*2-1,'-', 'Color',[0 0 0]+.4,'LineWidth',2.7);

for iT = 1:length(RespMat_logPR_cum{iN}(:,1))
    if RespMat{iN}(iT,2) == 1
        plot(iT, RespMat{iN}(iT,2), 'wo', 'MarkerFaceColor', cmap_c(1,:), 'MarkerSize', 9); 
    else
        plot(iT, RespMat{iN}(iT,2), 'wo', 'MarkerFaceColor', cmap_c(3,:), 'MarkerSize', 9); 
    end
%     if vectStable(iT) == 1
%         if RespMat{iN}(iT,2) == 1
%             plot(iT, RespMat{iN}(iT,2), 'wo', 'MarkerFaceColor', cmap_c(1,:), 'MarkerSize', 9); 
%         else
%             plot(iT, RespMat{iN}(iT,2), 'wo', 'MarkerFaceColor', cmap_c(3,:), 'MarkerSize', 9); 
%         end
%     else
%         if RespMat{iN}(iT,2) == 1
%             plot(iT, RespMat{iN}(iT,2)-0.1, 'o', 'MarkerFaceColor','w','Color', cmap_c(1,:), 'MarkerSize', 9,'linewidth',2); 
%         else
%             plot(iT, RespMat{iN}(iT,2)+0.1, 'o', 'MarkerFaceColor','w','Color', cmap_c(3,:), 'MarkerSize', 9,'linewidth',2); hold on;
%         end
%     end
end

yTickName{1}='ccw';yTickName{2}='cw';
set(gca,'yDir','normal');
 ylabel('Percept'); 
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [], 'YTick', [-1 1],'YTickLabel',yTickName, 'FontSize', Common_fontsize,'xcolor','w')
xlim([-1 102]); ylim([-1.3 1.3]); 





%% Change here to select the example trials
xtemp = 73;
iT_trial_analysis = [xtemp xtemp+1 xtemp+2]; 
hold on; 
for iT = 1:length(iT_trial_analysis)
    plot(iT_trial_analysis(iT), 1.3,'kv','color',[0 0 0]+ iT*0.2,'markerfacecolor',[0 0 0]+ iT*0.2, 'MarkerSize', 6); 
end
text(iT_trial_analysis(1), 1.56, {'Example trial #1~3'}, 'FontSize', Common_fontsize-2);







% LogPR evidence plot 
SP = subplot('position', [0.150    0.81    0.7    0.05]); cla; hold on;
yyaxis left

for iCase = 1:length(matRun(:,2))
    V_x = [matRun(iCase,4)-0.5 matRun(iCase,4):matRun(iCase,3) matRun(iCase,3)+0.5]; 
    V_y = [0 abs(RespMat_logPR_cum{iN}(matRun(iCase,4):matRun(iCase,3),end))' 0]; 
    V = [V_x ; V_y]'; 
    f = 1:length(V); 
    if matRun(iCase,1) == 1
        patch('Faces',f,'Vertices',V,'FaceColor',cmap_c(1,:),'edgecolor','none','FaceAlpha',0.8);
    else
        patch('Faces',f,'Vertices',V,'FaceColor',cmap_c(3,:),'edgecolor','none','FaceAlpha',0.8);
    end
end
xlim([-1 102]); ylim([-1.5 25]);
ylabel({'|logPR| ','& RT '});
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [0:20:100], 'XAxisLocation', 'bottom', 'YTick', [0 20], 'FontSize', Common_fontsize,'ycolor','k');

yyaxis right
hold on;
plot(find(vectStable==1), RT_model_plot(vectStable==1),'ko','markerfacecolor','k','markersize',4,'markeredgecolor',[0 0 0]+0.5);
plot(find(vectStable~=1), RT_model_plot(vectStable~=1),'ko','markerfacecolor','k','markersize',4,'markeredgecolor',[0 0 0]+0.5);
set(SP, 'YTick', [],'ycolor','k', 'FontSize', Common_fontsize);
ylim([-1.5 24]);
xlabel('Time (trial)'); 

% SP = subplot('position', [0.150    0.75    0.7    0.05]); cla; 
% yyaxis left
% hold on;
% for iT = 1:length(RespMat_logPR_cum{iN}(:,1))
%     if vectStable(iT) == 1
%         if RespMat{iN}(iT,2) == 1
%             bar(iT, abs(RespMat_logPR_cum{iN}(iT,end)),'FaceColor',cmap_c(3,:), 'EdgeColor',cmap_c(3,:),'Barwidth',0.8,'FaceAlpha',0.8); 
%         else
%             bar(iT, abs(RespMat_logPR_cum{iN}(iT,end)),'FaceColor',cmap_c(1,:), 'EdgeColor',cmap_c(1,:),'Barwidth',0.8,'FaceAlpha',0.8);
%         end
%     else
%         if RespMat{iN}(iT,2) == 1
%             bar(iT, abs(RespMat_logPR_cum{iN}(iT,end)),'FaceColor','w', 'EdgeColor',cmap_c(3,:),'Barwidth',0.8,'FaceAlpha',0.8); 
%         else
%             bar(iT, abs(RespMat_logPR_cum{iN}(iT,end)),'FaceColor','w', 'EdgeColor',cmap_c(1,:),'Barwidth',0.8,'FaceAlpha',0.8);
%         end
%     end
% end
% xlim([-1 102]); ylim([-1.5 25]);



SP = subplot('position', [0.880    0.81    0.03    0.05]); cla; hold on; 
line([3 4], [mean(RT_model_plot(vectStable==1)) mean(RT_model_plot(vectStable~=1))],'color',[0 0 0]+0.5); 
errorbar(3, mean(RT_model_plot(vectStable==1)), std(RT_model_plot(vectStable==1)),'ko','markerfacecolor','k','markeredgecolor',[0 0 0]+0.5,'color',[0 0 0]+0.5); 
errorbar(4, mean(RT_model_plot(vectStable~=1)), std(RT_model_plot(vectStable~=1)),'ko','markerfacecolor','w','markeredgecolor',[0 0 0]+0.5,'color',[0 0 0]+0.5); 
xlim([2.3 4.6]); ylim([-1.5 24]);
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [], 'YTick', [], 'FontSize', Common_fontsize,'ycolor','w')

SP = subplot('position', [0.920    0.81    0.03    0.05]); cla; hold on; 
line([1.0 2.0], [mean(RespMat_logPR_cum{iN}((vectStable==1) & (RespMat{iN}(:,2)==1)',end)) mean(RespMat_logPR_cum{iN}((vectStable~=1) & (RespMat{iN}(:,2)==1)',end))],'color',cmap_c(1,:)); 
errorbar(1.0, mean(RespMat_logPR_cum{iN}((vectStable==1) & (RespMat{iN}(:,2)==1)',end)), std(RespMat_logPR_cum{iN}((vectStable==1) & (RespMat{iN}(:,2)==1)',end)),'ko','markerfacecolor',cmap_c(1,:),'markeredgecolor',cmap_c(1,:),'color',cmap_c(3,:)); hold on; 
errorbar(2.0, mean(RespMat_logPR_cum{iN}((vectStable~=1) & (RespMat{iN}(:,2)==1)',end)), std(RespMat_logPR_cum{iN}((vectStable~=1) & (RespMat{iN}(:,2)==1)',end)),'ko','markerfacecolor','w','markeredgecolor',cmap_c(1,:),'color',cmap_c(1,:)); 
line([1.2 2.2], [-mean(RespMat_logPR_cum{iN}((vectStable==1) & (RespMat{iN}(:,2)~=1)',end)) -mean(RespMat_logPR_cum{iN}((vectStable~=1) & (RespMat{iN}(:,2)~=1)',end))],'color',cmap_c(3,:)); 
errorbar(1.2, -mean(RespMat_logPR_cum{iN}((vectStable==1) & (RespMat{iN}(:,2)~=1)',end)), std(RespMat_logPR_cum{iN}((vectStable==1) & (RespMat{iN}(:,2)~=1)',end)),'ko','markerfacecolor',cmap_c(3,:),'markeredgecolor',cmap_c(3,:),'color',cmap_c(3,:)); hold on; 
errorbar(2.2, -mean(RespMat_logPR_cum{iN}((vectStable~=1) & (RespMat{iN}(:,2)~=1)',end)), std(RespMat_logPR_cum{iN}((vectStable~=1) & (RespMat{iN}(:,2)~=1)',end)),'ko','markerfacecolor','w','markeredgecolor',cmap_c(3,:),'color',cmap_c(3,:)); 
xlim([0.6 2.6]); ylim([-1.5 24]);
% set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [1.1 2.1], 'XTicklabel',{'S', 'US'}, 'YTick', [], 'FontSize', 12,'xcolor','w','ycolor','w')
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [], 'YTick', [], 'FontSize', Common_fontsize,'xcolor','k','ycolor','w')







% likelihood 
SP = subplot('position', [0.15-0.1    0.72    0.09    0.06]); cla; 
plot(LKHDparam.axisVal, 0.2*0.8*LKHDparam.pdf_cw,'color',cmap_c(1,:),'linewidth',3); hold on; 
plot(LKHDparam.axisVal, 0.8*0.1*LKHDparam.pdf_ccw,'color',cmap_c(3,:),'linewidth',3);
line([0 0],[0 0.034],'linestyle','--','color','k'); 
xlim([min(LKHDparam.axisVal) max(LKHDparam.axisVal)]); ylim([-0.001 0.035]); 
xlabel('measurement (m)'); 
% title('Likelihood, p(m|\Omega)'); 
set(SP, 'box', 'off', 'TickDir', 'out','XTick',-5:5:5,'YTick',[],'FontSize', Common_fontsize,'ycolor','w');


% Prior (1st trial)
SP = subplot('position', [0.15-0.1    0.57    0.09    0.06]); cla; 
bar(1, result.matPrior(iT_trial_analysis(1)), 'facecolor', cmap_c(1,:) ,'edgecolor','w','barwidth',1); hold on; 
bar(-1, 1-result.matPrior(iT_trial_analysis(1)), 'facecolor', cmap_c(3,:) ,'edgecolor','w','barwidth',1);
line([-4 4], [0 0]+0.5,'color','k','linestyle','-.'); 
xlim([-2.5 2.5]); ylim([0 1]); %ylabel('probability'); 
text(1-0.7, 0.9, num2str(round(result.matPrior(iT_trial_analysis(1))*100)/100),'FontSize', Common_fontsize-2);
text(-1-0.7, 0.9, num2str(1 - round(result.matPrior(iT_trial_analysis(1))*100)/100),'FontSize', Common_fontsize-2);
xlabel('world state (\Omega)'); 
set(SP, 'box', 'off', 'TickDir', 'out','XTick',[-1 1],'XTicklabel',{'ccw','cw'},'YTick',[],'ycolor','w','FontSize', Common_fontsize);


% Posterior (1st trial)
SP = subplot('position', [0.15-0.1    0.45    0.09    0.06]); cla; 
% plot(LKHDparam.axisVal, result.matPrior(iT_trial_analysis(1))*LKHDparam.pdf_cw,'color',cmap_c(1,:),'linewidth',2.5); hold on; 
% plot(LKHDparam.axisVal, (1-result.matPrior(iT_trial_analysis(1)))*LKHDparam.pdf_ccw,'color',cmap_c(3,:),'linewidth',2.5);
plot(LKHDparam.axisVal, result.matPrior(iT_trial_analysis(1))*LKHDparam.pdf_cw + (1-result.matPrior(iT_trial_analysis(1)))*LKHDparam.pdf_ccw,'color','k','linewidth',2.5);
line([0 0],[0 0.034],'linestyle','--','color','k'); 
xlim([min(LKHDparam.axisVal) max(LKHDparam.axisVal)]); ylim([-0.001 0.03]) 
set(SP, 'box', 'off', 'TickDir', 'out','XTick',-5:5:5,'YTick',[],'FontSize', Common_fontsize,'ycolor','w');
% xlabel('measurement (m)');


SP = subplot('position', [0.15-0.1    0.22    0.07    0.2]); cla; 
% yyaxis left; 
hold on; 
for iSample = 1:POparam.sampleRep
    plot(iSample, RespMat_allsample{iN}(iT_trial_analysis(1),iSample), 'ko','markerfacecolor',[0 0 0]+0.7, 'markeredgecolor',[0 0 0]+0.7,'markersize',6);
%     if RespMat_whichDist{iN}(iT_trial_analysis(1),iSample) == 1
%         plot(iSample, RespMat_allsample{iN}(iT_trial_analysis(1),iSample), 'ko','markerfacecolor',[0 0 0]+0.7, 'markeredgecolor',cmap_c(1,:),'markersize',6);
%     else
%         plot(iSample, RespMat_allsample{iN}(iT_trial_analysis(1),iSample), 'ko','markerfacecolor',[0 0 0]+0.7, 'markeredgecolor',cmap_c(3,:),'markersize',6);
%     end
end
line([0 POparam.sampleSet(iN)+1],[0 0],'linestyle','--','color','k'); 
% plot(1:iSample, vectOutcome_A(1:iSample),'ko','markerfacecolor',[0 0 0]+0.3, 'markeredgecolor','w','markersize',8);
xlim([0 POparam.sampleSet(iN)+1]); ylim([-6 6]); 
ylabel('m of sample');
xlabel('within-trial time (samples)'); 
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[0:5:20],'Ytick',[-5 0 5],'ycolor','k');
view(90,90); 



SP = subplot('position', [0.15-0.1    0.11    0.09    0.06]); cla; hold on; 
M_range = -7:1:7;
[cwMat,~] = hist(RespMat_allsample{iN}(iT_trial_analysis(1),RespMat_whichDist{iN}(iT_trial_analysis(1),:)==1), M_range); 
[ccwMat,~] = hist(RespMat_allsample{iN}(iT_trial_analysis(1),RespMat_whichDist{iN}(iT_trial_analysis(1),:)~=1), M_range); 
h = bar(M_range, [cwMat; ccwMat]','stacked'); 
set(h(1), 'FaceColor',cmap_c(2,:),'BarWidth',0.8,'EdgeColor','none'); 
set(h(2), 'FaceColor',cmap_c(2,:),'BarWidth',0.8,'EdgeColor','none'); 
line([0 0],[0 6],'linestyle','--','color','k'); 
plot(mean(RespMat_allsample{iN}(iT_trial_analysis(1),:)),0.5,'kv','markerfacecolor',cmap_c(1,:))
ylim([0 6]); xlim([-7 7]);
xlabel(''); 
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[-5:5:5],'Ytick',[0 2 4 6],'ycolor','k');



SP = subplot('position', [0.25-0.1-0.01    0.22    0.07    0.2]); cla; hold on; 
Vspace = [0 -POparam.DecisionBound; 20.9 -POparam.DecisionBound; 20.9 -16; 0 -16];
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(3,:),'edgecolor','w','facealpha',0.6); hold on; 
Vspace = [0 POparam.DecisionBound; 20.9 POparam.DecisionBound; 20.9 16; 0 16];
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(1,:),'edgecolor','w','facealpha',0.6); hold on; 
Vspace = [0 -POparam.DecisionBound; 20.9 -POparam.DecisionBound; 20.9 0; 0 0];
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(3,:),'edgecolor','w','facealpha',0.2); hold on; 
Vspace = [0 POparam.DecisionBound; 20.9 POparam.DecisionBound; 20.9 0; 0 0];
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(1,:),'edgecolor','w','facealpha',0.2); hold on; 
ylabel('logPR'); 
stairs((0:length(RespMat_logPR_cum{iN}(iT_trial_analysis(1),:))), [0 RespMat_logPR_cum{iN}(iT_trial_analysis(1),:)],'color',[0 0 0]+0.7,'linewidth',2.5); hold on; %plot(1:length(SumMat), SumMat,'ko','markerfacecolor','w','markeredgecolor','k');
stairs((0:RT_model(iT_trial_analysis(1))), [0 RespMat_logPR_cum{iN}(iT_trial_analysis(1),1:RT_model(iT_trial_analysis(1)))],'color',[0 0 0]+0.1,'linewidth',2.5); hold on; %plot(1:length(SumMat), SumMat,'ko','markerfacecolor','w','markeredgecolor','k');
for iSample = 1:POparam.sampleRep
    if RespMat_logPR{iN}(iT_trial_analysis(1),iSample)>0
        line([iSample, iSample], [RespMat_logPR_cum{iN}(iT_trial_analysis(1),iSample) (RespMat_logPR_cum{iN}(iT_trial_analysis(1),iSample)-RespMat_logPR{iN}(iT_trial_analysis(1),iSample))], 'color',cmap_c(1,:),'linewidth',3); 
    else
        line([iSample, iSample], [RespMat_logPR_cum{iN}(iT_trial_analysis(1),iSample) (RespMat_logPR_cum{iN}(iT_trial_analysis(1),iSample)-RespMat_logPR{iN}(iT_trial_analysis(1),iSample))], 'color',cmap_c(3,:),'linewidth',3); 
    end
end
if RespMat{iN}(iT_trial_analysis(1),2) == 1
    plot(RT_model(iT_trial_analysis(1)), RespMat_logPR_cum{iN}(iT_trial_analysis(1),RT_model(iT_trial_analysis(1))),'k>','markerfacecolor',cmap_c(1,:),'markersize',7);
else
    plot(RT_model(iT_trial_analysis(1)), RespMat_logPR_cum{iN}(iT_trial_analysis(1),RT_model(iT_trial_analysis(1))),'k<','markerfacecolor',cmap_c(3,:),'markersize',7);
end
ylim([-16 16]); xlim([0 POparam.sampleSet(iN)+1]);
line([0 POparam.sampleSet(iN)+1],[0 0],'linestyle','--','color','k'); 
view(90,90); 
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[],'Ytick',[-10 0 10],'xcolor','w');



%% Trial 1 --> 2 (update)

SP = subplot('position', [0.15-0.1+0.2+0.01    0.57    0.045    0.032]); cla; hold on; 
plot(PRparam.vectPi, matPrior(iT_trial_analysis(1),:)/max(matPrior(iT_trial_analysis(1),:)),'color','k','linewidth',1);
if RespMat{iN}(iT_trial_analysis(1),2) == 1
    plot(PRparam.vectPi, matLKHD(iT_trial_analysis(1),:)/max(matLKHD(iT_trial_analysis(1),:)),'color',cmap_c(1,:),'linewidth',1); 
else
    plot(PRparam.vectPi, matLKHD(iT_trial_analysis(1),:)/max(matLKHD(iT_trial_analysis(1),:)),'color',cmap_c(3,:),'linewidth',1); 
end
plot(PRparam.vectPi, matPOST(iT_trial_analysis(1),:)/max(matPOST(iT_trial_analysis(1),:)),'color',[0 0.6 0],'linewidth',1);
ylim([-0.02 1.1])
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize-2,'Xtick',[0 1],'Ytick',[],'ycolor','k');
yyaxis right
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize-2,'Xtick',[0 1],'Ytick',[],'ycolor','k');
xlabel('\omega');

staticBase = (1-PRparam.alpha)/sum(PRparam.vectPi); 

SP = subplot('position', [0.15-0.1+0.2+0.07    0.57    0.045    0.032]); cla; hold on; 
plot(PRparam.vectPi, matPOST(iT_trial_analysis(1),:)*PRparam.alpha + staticBase,'color',[0 0.6 0],'linewidth',1); 
line([0 1],[0 0]+staticBase,'color','k','linestyle','-.','linewidth',1);
tempPost = matPOST(iT_trial_analysis(1),:)*PRparam.alpha + staticBase; 
tempX = [PRparam.vectPi(1:find(PRparam.vectPi == 0.5)); tempPost(1:find(PRparam.vectPi == 0.5))]';
Vspace = [tempX; 0.5 0; 0 0; 0 tempPost(1)]; 
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(3,:),'edgecolor','none','facealpha',0.5); hold on; 
tempX = [PRparam.vectPi(find(PRparam.vectPi == 0.5):end); tempPost(find(PRparam.vectPi == 0.5):end)]';
Vspace = [tempX; 1 0; 0.5 0 ; 0.5 tempPost((PRparam.vectPi == 0.5));]; 
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(1,:),'edgecolor','none','facealpha',0.5); 
xlabel('\omega');
xlim([-0.02 1.02]); ylim([0 max((matPOST(iT_trial_analysis(1),:)*PRparam.alpha + staticBase))*1.1]); 
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize-2,'Xtick',[0 1],'Ytick',[],'ycolor','k');
yyaxis right
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize-2,'Xtick',[0 1],'Ytick',[],'ycolor','k');

% hold on; 
% 
% staticBase = (1-PRparam.alpha)/sum(PRparam.vectPi); 
% tempPost = matPOST(iT_trial_analysis(1),:)*PRparam.alpha + staticBase; 
% tempX = [PRparam.vectPi(1:find(PRparam.vectPi == 0.5)); tempPost(1:find(PRparam.vectPi == 0.5))]';
% Vspace = [tempX; 0.5 0; 0 0; 0 tempPost(1)]; 
% f = 1:length(Vspace);
% patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(3,:),'edgecolor','w','facealpha',0.4); hold on; 
% tempX = [PRparam.vectPi(find(PRparam.vectPi == 0.5):end); tempPost(find(PRparam.vectPi == 0.5):end)]';
% Vspace = [tempX; 1 0; 0.5 0 ; 0.5 tempPost((PRparam.vectPi == 0.5));]; 
% f = 1:length(Vspace);
% patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(1,:),'edgecolor','w','facealpha',0.4); 
% 
% plot(PRparam.vectPi, matPrior(iT_trial_analysis(1),:)*PRparam.alpha + staticBase,'color','k','linewidth',1.5); 
% if RespMat{iN}(iT_trial_analysis(1),2) == 1
%     plot(PRparam.vectPi, matLKHD(iT_trial_analysis(1),:)/sum(matLKHD(iT_trial_analysis(1),:))*PRparam.alpha + staticBase,'color',cmap_c(1,:),'linewidth',2); 
% else
%     plot(PRparam.vectPi, matLKHD(iT_trial_analysis(1),:)/sum(matLKHD(iT_trial_analysis(1),:))*PRparam.alpha + staticBase,'color',cmap_c(3,:),'linewidth',2); 
% end
% plot(PRparam.vectPi, matPOST(iT_trial_analysis(1),:)*PRparam.alpha + staticBase,'color',[0 0.6 0],'linewidth',2); 
% line([0 1],[0 0]+staticBase,'color','k','linestyle','-.','linewidth',1.5); 
% ValueAll = [(matPrior(iT_trial_analysis(1),:)*PRparam.alpha + staticBase) (matPOST(iT_trial_analysis(1),:)*PRparam.alpha + staticBase) (matLKHD(iT_trial_analysis(1),:)/sum(matLKHD(iT_trial_analysis(1),:))*PRparam.alpha + staticBase)];
% xlim([-0.02 1.02]); ylim([0 max(ValueAll)*1.01]); xlabel('\omega');
% set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[0 0.5 1],'Ytick',[],'ycolor','k');
% yyaxis right
% set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[0 0.5 1],'Ytick',[],'ycolor','k');
%% Trial 2


% Prior (1st trial)
SP = subplot('position', [0.48+0.02-0.1+0.01    0.57    0.09    0.06]); cla; 
bar(1, result.matPrior(iT_trial_analysis(2)), 'facecolor', cmap_c(1,:) ,'edgecolor','w','barwidth',1); hold on; 
bar(-1, 1-result.matPrior(iT_trial_analysis(2)), 'facecolor', cmap_c(3,:) ,'edgecolor','w','barwidth',1);
line([-4 4], [0 0]+0.5,'color','k','linestyle','-.'); 
xlim([-2.5 2.5]); ylim([0 1]); %ylabel('probability'); 
text(1-0.7, 0.9, num2str(round(result.matPrior(iT_trial_analysis(2))*100)/100),'FontSize', Common_fontsize-2);
text(-1-0.7, 0.9, num2str(1 - round(result.matPrior(iT_trial_analysis(2))*100)/100),'FontSize', Common_fontsize-2);
xlabel('world state (\Omega)'); 
set(SP, 'box', 'off', 'TickDir', 'out','XTick',[-1 1],'XTicklabel',{'ccw','cw'},'YTick',[],'ycolor','w','FontSize', Common_fontsize);


% Posterior (1st trial)
SP = subplot('position', [0.48+0.02-0.1+0.01    0.45    0.09    0.06]); cla; 
% plot(LKHDparam.axisVal, result.matPrior(iT_trial_analysis(2))*LKHDparam.pdf_cw,'color',cmap_c(1,:),'linewidth',2.5); hold on; 
% plot(LKHDparam.axisVal, (1-result.matPrior(iT_trial_analysis(2)))*LKHDparam.pdf_ccw,'color',cmap_c(3,:),'linewidth',2.5);
plot(LKHDparam.axisVal, result.matPrior(iT_trial_analysis(2))*LKHDparam.pdf_cw + (1-result.matPrior(iT_trial_analysis(2)))*LKHDparam.pdf_ccw,'color','k','linewidth',2.5);
line([0 0],[0 0.034],'linestyle','--','color','k'); 
xlim([min(LKHDparam.axisVal) max(LKHDparam.axisVal)]); ylim([-0.001 0.03]) 
set(SP, 'box', 'off', 'TickDir', 'out','XTick',-5:5:5,'YTick',[],'FontSize', Common_fontsize,'ycolor','w');
% xlabel('measurement (m)');


SP = subplot('position', [0.48+0.02-0.1+0.01    0.22    0.07    0.2]); cla; 
% yyaxis left; 
hold on; 
for iSample = 1:POparam.sampleRep
    plot(iSample, RespMat_allsample{iN}(iT_trial_analysis(2),iSample), 'ko','markerfacecolor',[0 0 0]+0.7, 'markeredgecolor',[0 0 0]+0.7,'markersize',6);
%     if RespMat_whichDist{iN}(iT_trial_analysis(2),iSample) == 1
%         plot(iSample, RespMat_allsample{iN}(iT_trial_analysis(2),iSample), 'ko','markerfacecolor',[0 0 0]+0.7, 'markeredgecolor',cmap_c(1,:),'markersize',6);
%     else
%         plot(iSample, RespMat_allsample{iN}(iT_trial_analysis(2),iSample), 'ko','markerfacecolor',[0 0 0]+0.7, 'markeredgecolor',cmap_c(3,:),'markersize',6);
%     end
end
line([0 POparam.sampleSet(iN)+1],[0 0],'linestyle','--','color','k'); 
% plot(1:iSample, vectOutcome_A(1:iSample),'ko','markerfacecolor',[0 0 0]+0.3, 'markeredgecolor','w','markersize',8);
xlim([0 POparam.sampleSet(iN)+1]); ylim([-6 6]); 
ylabel('m of sample');
xlabel('within-trial time (samples)'); 
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[0:5:20],'Ytick',[-5 0 5],'ycolor','k');
view(90,90); 



SP = subplot('position', [0.48+0.02-0.1+0.01    0.11    0.09    0.06]); cla; hold on; 
M_range = -7:1:7;
[cwMat,~] = hist(RespMat_allsample{iN}(iT_trial_analysis(2),RespMat_whichDist{iN}(iT_trial_analysis(2),:)==1), M_range); 
[ccwMat,~] = hist(RespMat_allsample{iN}(iT_trial_analysis(2),RespMat_whichDist{iN}(iT_trial_analysis(2),:)~=1), M_range); 
h = bar(M_range, [cwMat; ccwMat]','stacked'); 
set(h(1), 'FaceColor',cmap_c(2,:),'BarWidth',0.8,'EdgeColor','none'); 
set(h(2), 'FaceColor',cmap_c(2,:),'BarWidth',0.8,'EdgeColor','none'); 
line([0 0],[0 6],'linestyle','--','color','k'); 
plot(mean(RespMat_allsample{iN}(iT_trial_analysis(2),:)),0.5,'kv','markerfacecolor',cmap_c(3,:))
ylim([0 6]); xlim([-7 7]);
xlabel(''); 
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[-5:5:5],'Ytick',[0 2 4 6],'ycolor','k');



SP = subplot('position', [0.27+0.35-0.1-0.02    0.22    0.07    0.2]); cla; hold on; 
Vspace = [0 -POparam.DecisionBound; 20.9 -POparam.DecisionBound; 20.9 -16; 0 -16];
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(3,:),'edgecolor','w','facealpha',0.6); hold on; 
Vspace = [0 POparam.DecisionBound; 20.9 POparam.DecisionBound; 20.9 16; 0 16];
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(1,:),'edgecolor','w','facealpha',0.6); hold on; 
Vspace = [0 -POparam.DecisionBound; 20.9 -POparam.DecisionBound; 20.9 0; 0 0];
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(3,:),'edgecolor','w','facealpha',0.2); hold on; 
Vspace = [0 POparam.DecisionBound; 20.9 POparam.DecisionBound; 20.9 0; 0 0];
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(1,:),'edgecolor','w','facealpha',0.2); hold on; 
ylabel('logPR'); 
stairs((0:length(RespMat_logPR_cum{iN}(iT_trial_analysis(2),:))), [0 RespMat_logPR_cum{iN}(iT_trial_analysis(2),:)],'color',[0 0 0]+0.7,'linewidth',2.5); hold on; %plot(1:length(SumMat), SumMat,'ko','markerfacecolor','w','markeredgecolor','k');
stairs((0:RT_model(iT_trial_analysis(2))), [0 RespMat_logPR_cum{iN}(iT_trial_analysis(2),1:RT_model(iT_trial_analysis(2)))],'color',[0 0 0]+0.1,'linewidth',2.5); hold on; %plot(1:length(SumMat), SumMat,'ko','markerfacecolor','w','markeredgecolor','k');
for iSample = 1:POparam.sampleRep
    if RespMat_logPR{iN}(iT_trial_analysis(2),iSample)>0
        line([iSample, iSample], [RespMat_logPR_cum{iN}(iT_trial_analysis(2),iSample) (RespMat_logPR_cum{iN}(iT_trial_analysis(2),iSample)-RespMat_logPR{iN}(iT_trial_analysis(2),iSample))], 'color',cmap_c(1,:),'linewidth',3); 
    else
        line([iSample, iSample], [RespMat_logPR_cum{iN}(iT_trial_analysis(2),iSample) (RespMat_logPR_cum{iN}(iT_trial_analysis(2),iSample)-RespMat_logPR{iN}(iT_trial_analysis(2),iSample))], 'color',cmap_c(3,:),'linewidth',3); 
    end
end
if RespMat{iN}(iT_trial_analysis(2),2) == 1
    plot(RT_model(iT_trial_analysis(2)), RespMat_logPR_cum{iN}(iT_trial_analysis(2),RT_model(iT_trial_analysis(2))),'k>','markerfacecolor',cmap_c(1,:),'markersize',7);
else
    plot(RT_model(iT_trial_analysis(2)), RespMat_logPR_cum{iN}(iT_trial_analysis(2),RT_model(iT_trial_analysis(2))),'k<','markerfacecolor',cmap_c(3,:),'markersize',7);
end
ylim([-16 16]); xlim([0 POparam.sampleSet(iN)+1]);
line([0 POparam.sampleSet(iN)+1],[0 0],'linestyle','--','color','k'); 
view(90,90); 
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[],'Ytick',[-10 0 10],'xcolor','w');



%% Trial 2 --> 3 (update)

SP = subplot('position', [0.15-0.1+0.2+0.01+0.35    0.57    0.045    0.032]); cla; hold on; 
plot(PRparam.vectPi, matPrior(iT_trial_analysis(2),:)/max(matPrior(iT_trial_analysis(2),:)),'color','k','linewidth',1);
if RespMat{iN}(iT_trial_analysis(2),2) == 1
    plot(PRparam.vectPi, matLKHD(iT_trial_analysis(2),:)/max(matLKHD(iT_trial_analysis(2),:)),'color',cmap_c(1,:),'linewidth',1); 
else
    plot(PRparam.vectPi, matLKHD(iT_trial_analysis(2),:)/max(matLKHD(iT_trial_analysis(2),:)),'color',cmap_c(3,:),'linewidth',1); 
end
plot(PRparam.vectPi, matPOST(iT_trial_analysis(2),:)/max(matPOST(iT_trial_analysis(2),:)),'color',[0 0.6 0],'linewidth',1);
ylim([-0.02 1.1])
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize-2,'Xtick',[0 1],'Ytick',[],'ycolor','k');
yyaxis right
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize-2,'Xtick',[0 1],'Ytick',[],'ycolor','k');

SP = subplot('position', [0.15-0.1+0.2+0.07+0.35    0.57    0.045    0.032]); cla; hold on; 
plot(PRparam.vectPi, matPOST(iT_trial_analysis(2),:)*PRparam.alpha + staticBase,'color',[0 0.6 0],'linewidth',1); 
line([0 1],[0 0]+staticBase,'color','k','linestyle','-.','linewidth',1);
staticBase = (1-PRparam.alpha)/sum(PRparam.vectPi); 
tempPost = matPOST(iT_trial_analysis(2),:)*PRparam.alpha + staticBase; 
tempX = [PRparam.vectPi(1:find(PRparam.vectPi == 0.5)); tempPost(1:find(PRparam.vectPi == 0.5))]';
Vspace = [tempX; 0.5 0; 0 0; 0 tempPost(1)]; 
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(3,:),'edgecolor','none','facealpha',0.5); hold on; 
tempX = [PRparam.vectPi(find(PRparam.vectPi == 0.5):end); tempPost(find(PRparam.vectPi == 0.5):end)]';
Vspace = [tempX; 1 0; 0.5 0 ; 0.5 tempPost((PRparam.vectPi == 0.5));]; 
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(1,:),'edgecolor','none','facealpha',0.5); 
xlim([-0.02 1.02]); ylim([0 max((matPOST(iT_trial_analysis(2),:)*PRparam.alpha + staticBase))*1.1]); xlabel('\omega');
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize-2,'Xtick',[0 1],'Ytick',[],'ycolor','k');
yyaxis right
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize-2,'Xtick',[0 1],'Ytick',[],'ycolor','k');



% % SP = subplot('position', [0.15-0.1+0.2+0.35+0.02    0.54    0.07    0.05]); cla; 
% SP = subplot('position', [0.15-0.1+0.2+0.35+0.01    0.57    0.045    0.032]); cla; hold on; 
% plot(PRparam.vectPi, matPrior(iT_trial_analysis(2),:)/max(matPrior(iT_trial_analysis(2),:)),'color','k','linewidth',1.5);
% if RespMat{iN}(iT_trial_analysis(2),2) == 1
%     plot(PRparam.vectPi, matLKHD(iT_trial_analysis(2),:)/max(matLKHD(iT_trial_analysis(2),:)),'color',cmap_c(1,:),'linewidth',2); 
% else
%     plot(PRparam.vectPi, matLKHD(iT_trial_analysis(2),:)/max(matLKHD(iT_trial_analysis(2),:)),'color',cmap_c(3,:),'linewidth',2); 
% end
% plot(PRparam.vectPi, matPOST(iT_trial_analysis(2),:)/max(matPOST(iT_trial_analysis(2),:)),'color',[0 0.6 0],'linewidth',2);
% ylim([-0.02 1.1])
% set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize-2,'Xtick',[0 1],'Ytick',[],'ycolor','k');
% yyaxis right
% set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize-2,'Xtick',[0 1],'Ytick',[],'ycolor','k');
% 
% SP = subplot('position', [0.15-0.1+0.2+0.35+0.07    0.57    0.045    0.032]); cla; hold on; 
% staticBase = (1-PRparam.alpha)/sum(PRparam.vectPi); 
% tempPost = matPOST(iT_trial_analysis(2),:)*PRparam.alpha + staticBase; 
% tempX = [PRparam.vectPi(1:find(PRparam.vectPi == 0.5)); tempPost(1:find(PRparam.vectPi == 0.5))]';
% Vspace = [tempX; 0.5 0; 0 0; 0 tempPost(1)]; 
% f = 1:length(Vspace);
% patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(3,:),'edgecolor','w','facealpha',0.4); hold on; 
% tempX = [PRparam.vectPi(find(PRparam.vectPi == 0.5):end); tempPost(find(PRparam.vectPi == 0.5):end)]';
% Vspace = [tempX; 1 0; 0.5 0 ; 0.5 tempPost((PRparam.vectPi == 0.5));]; 
% f = 1:length(Vspace);
% patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(1,:),'edgecolor','w','facealpha',0.4); 
% plot(PRparam.vectPi, matPOST(iT_trial_analysis(2),:)*PRparam.alpha + staticBase,'color',[0 0.6 0],'linewidth',2); 
% line([0 1],[0 0]+staticBase,'color','k','linestyle','-.','linewidth',1.5);
% xlim([-0.02 1.02]); ylim([0 max((matPOST(iT_trial_analysis(2),:)*PRparam.alpha + staticBase))*1.1]); xlabel('\omega');
% set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize-2,'Xtick',[0 1],'Ytick',[],'ycolor','k');
% yyaxis right
% set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize-2,'Xtick',[0 1],'Ytick',[],'ycolor','k');
% 
% % hold on; 
% % 
% % staticBase = (1-PRparam.alpha)/sum(PRparam.vectPi); 
% % tempPost = matPOST(iT_trial_analysis(2),:)*PRparam.alpha + staticBase; 
% % tempX = [PRparam.vectPi(1:find(PRparam.vectPi == 0.5)); tempPost(1:find(PRparam.vectPi == 0.5))]';
% % Vspace = [tempX; 0.5 0; 0 0; 0 tempPost(1)]; 
% % f = 1:length(Vspace);
% % patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(3,:),'edgecolor','w','facealpha',0.4); hold on; 
% % tempX = [PRparam.vectPi(find(PRparam.vectPi == 0.5):end); tempPost(find(PRparam.vectPi == 0.5):end)]';
% % Vspace = [tempX; 1 0; 0.5 0 ; 0.5 tempPost((PRparam.vectPi == 0.5));]; 
% % f = 1:length(Vspace);
% % patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(1,:),'edgecolor','w','facealpha',0.4); 
% % 
% % plot(PRparam.vectPi, matPrior(iT_trial_analysis(2),:)*PRparam.alpha + staticBase,'color','k','linewidth',1.5); 
% % if RespMat{iN}(iT_trial_analysis(2),2) == 1
% %     plot(PRparam.vectPi, matLKHD(iT_trial_analysis(2),:)/sum(matLKHD(iT_trial_analysis(2),:))*PRparam.alpha + staticBase,'color',cmap_c(1,:),'linewidth',2); 
% % else
% %     plot(PRparam.vectPi, matLKHD(iT_trial_analysis(2),:)/sum(matLKHD(iT_trial_analysis(2),:))*PRparam.alpha + staticBase,'color',cmap_c(3,:),'linewidth',2); 
% % end
% % plot(PRparam.vectPi, matPOST(iT_trial_analysis(2),:)*PRparam.alpha + staticBase,'color',[0 0.6 0],'linewidth',2); 
% % line([0 1],[0 0]+staticBase,'color','k','linestyle','-.','linewidth',1.5); 
% % ValueAll = [(matPrior(iT_trial_analysis(2),:)*PRparam.alpha + staticBase) (matPOST(iT_trial_analysis(2),:)*PRparam.alpha + staticBase) (matLKHD(iT_trial_analysis(2),:)/sum(matLKHD(iT_trial_analysis(2),:))*PRparam.alpha + staticBase)];
% % xlim([-0.02 1.02]); ylim([0 max(ValueAll)*1.01]); xlabel('\omega');
% % set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[0 0.5 1],'Ytick',[],'ycolor','k');
% % yyaxis right
% % set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[0 0.5 1],'Ytick',[],'ycolor','k');




%% Trial 3


% Prior (1st trial)
SP = subplot('position', [0.48+0.02+0.35+0.02-0.1    0.57    0.09    0.06]); cla; 
bar(1, result.matPrior(iT_trial_analysis(3)), 'facecolor', cmap_c(1,:) ,'edgecolor','w','barwidth',1); hold on; 
bar(-1, 1-result.matPrior(iT_trial_analysis(3)), 'facecolor', cmap_c(3,:) ,'edgecolor','w','barwidth',1);
line([-4 4], [0 0]+0.5,'color','k','linestyle','-.'); 
xlim([-2.5 2.5]); ylim([0 1]); %ylabel('probability'); 
text(1-0.4, 0.9, num2str(round(result.matPrior(iT_trial_analysis(3))*100)/100),'FontSize', Common_fontsize-2);
text(-1-0.4, 0.9, num2str(1 - round(result.matPrior(iT_trial_analysis(3))*100)/100),'FontSize', Common_fontsize-2);
xlabel('world state (\Omega)'); 
set(SP, 'box', 'off', 'TickDir', 'out','XTick',[-1 1],'XTicklabel',{'ccw','cw'},'YTick',[],'ycolor','w','FontSize', Common_fontsize);


% Posterior (1st trial)
SP = subplot('position', [0.48+0.02+0.35+0.02-0.1    0.45    0.09    0.06]); cla; 
% plot(LKHDparam.axisVal, result.matPrior(iT_trial_analysis(3))*LKHDparam.pdf_cw,'color',cmap_c(1,:),'linewidth',2.5); hold on; 
% plot(LKHDparam.axisVal, (1-result.matPrior(iT_trial_analysis(3)))*LKHDparam.pdf_ccw,'color',cmap_c(3,:),'linewidth',2.5);
plot(LKHDparam.axisVal, result.matPrior(iT_trial_analysis(3))*LKHDparam.pdf_cw + (1-result.matPrior(iT_trial_analysis(3)))*LKHDparam.pdf_ccw,'color','k','linewidth',2.5);
line([0 0],[0 0.034],'linestyle','--','color','k'); 
xlim([min(LKHDparam.axisVal) max(LKHDparam.axisVal)]); ylim([-0.001 0.03]) 
set(SP, 'box', 'off', 'TickDir', 'out','XTick',-5:5:5,'YTick',[],'FontSize', Common_fontsize,'ycolor','w');
% xlabel('measurement (m)');









SP = subplot('position', [0.48+0.02+0.35-0.1+0.02   0.22    0.07    0.2]); cla; 
% yyaxis left; 
hold on; 
for iSample = 1:POparam.sampleRep
    plot(iSample, RespMat_allsample{iN}(iT_trial_analysis(3),iSample), 'ko','markerfacecolor',[0 0 0]+0.7, 'markeredgecolor',[0 0 0]+0.7,'markersize',6);
%     if RespMat_whichDist{iN}(iT_trial_analysis(3),iSample) == 1
%         plot(iSample, RespMat_allsample{iN}(iT_trial_analysis(3),iSample), 'ko','markerfacecolor',[0 0 0]+0.7, 'markeredgecolor',cmap_c(1,:),'markersize',6);
%     else
%         plot(iSample, RespMat_allsample{iN}(iT_trial_analysis(3),iSample), 'ko','markerfacecolor',[0 0 0]+0.7, 'markeredgecolor',cmap_c(3,:),'markersize',6);
%     end
end
line([0 POparam.sampleSet(iN)+1],[0 0],'linestyle','--','color','k'); 
% plot(1:iSample, vectOutcome_A(1:iSample),'ko','markerfacecolor',[0 0 0]+0.3, 'markeredgecolor','w','markersize',8);
xlim([0 POparam.sampleSet(iN)+1]); ylim([-6 6]); 
ylabel('m of sample');
xlabel('within-trial time (samples)'); 
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[0:5:20],'Ytick',[-5 0 5],'ycolor','k');
view(90,90); 



SP = subplot('position', [0.48+0.04+0.35-0.1    0.11    0.09    0.06]); cla; hold on; 
M_range = -7:1:7;
[cwMat,~] = hist(RespMat_allsample{iN}(iT_trial_analysis(3),RespMat_whichDist{iN}(iT_trial_analysis(3),:)==1), M_range); 
[ccwMat,~] = hist(RespMat_allsample{iN}(iT_trial_analysis(3),RespMat_whichDist{iN}(iT_trial_analysis(3),:)~=1), M_range); 
h = bar(M_range, [cwMat; ccwMat]','stacked'); 
set(h(1), 'FaceColor',cmap_c(2,:),'BarWidth',0.8,'EdgeColor','none'); 
set(h(2), 'FaceColor',cmap_c(2,:),'BarWidth',0.8,'EdgeColor','none'); 
line([0 0],[0 6],'linestyle','--','color','k'); 
plot(mean(RespMat_allsample{iN}(iT_trial_analysis(3),:)),0.5,'kv','markerfacecolor',cmap_c(1,:))
ylim([0 6]); xlim([-7 7]);
xlabel(''); 
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[-5:5:5],'Ytick',[0 2 4 6],'ycolor','k');



SP = subplot('position', [0.27+0.35+0.35-0.11    0.22    0.07    0.2]); cla; hold on; 
Vspace = [0 -POparam.DecisionBound; 20.9 -POparam.DecisionBound; 20.9 -16; 0 -16];
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(3,:),'edgecolor','w','facealpha',0.6); hold on; 
Vspace = [0 POparam.DecisionBound; 20.9 POparam.DecisionBound; 20.9 16; 0 16];
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(1,:),'edgecolor','w','facealpha',0.6); hold on; 
Vspace = [0 -POparam.DecisionBound; 20.9 -POparam.DecisionBound; 20.9 0; 0 0];
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(3,:),'edgecolor','w','facealpha',0.2); hold on; 
Vspace = [0 POparam.DecisionBound; 20.9 POparam.DecisionBound; 20.9 0; 0 0];
f = 1:length(Vspace);
patch('Faces',f,'Vertices',Vspace,'FaceColor',cmap_c(1,:),'edgecolor','w','facealpha',0.2); hold on; 
ylabel('logPR'); 

stairs((0:length(RespMat_logPR_cum{iN}(iT_trial_analysis(3),:))), [0 RespMat_logPR_cum{iN}(iT_trial_analysis(3),:)],'color',[0 0 0]+0.7,'linewidth',2.5); hold on; %plot(1:length(SumMat), SumMat,'ko','markerfacecolor','w','markeredgecolor','k');
stairs((0:RT_model(iT_trial_analysis(3))), [0 RespMat_logPR_cum{iN}(iT_trial_analysis(3),1:RT_model(iT_trial_analysis(3)))],'color',[0 0 0]+0.1,'linewidth',2.5); hold on; %plot(1:length(SumMat), SumMat,'ko','markerfacecolor','w','markeredgecolor','k');
for iSample = 1:POparam.sampleRep
    if RespMat_logPR{iN}(iT_trial_analysis(3),iSample)>0
        line([iSample, iSample], [RespMat_logPR_cum{iN}(iT_trial_analysis(3),iSample) (RespMat_logPR_cum{iN}(iT_trial_analysis(3),iSample)-RespMat_logPR{iN}(iT_trial_analysis(3),iSample))], 'color',cmap_c(1,:),'linewidth',3); 
    else
        line([iSample, iSample], [RespMat_logPR_cum{iN}(iT_trial_analysis(3),iSample) (RespMat_logPR_cum{iN}(iT_trial_analysis(3),iSample)-RespMat_logPR{iN}(iT_trial_analysis(3),iSample))], 'color',cmap_c(3,:),'linewidth',3); 
    end
end
if RespMat{iN}(iT_trial_analysis(3),2) == 1
    plot(RT_model(iT_trial_analysis(3)), RespMat_logPR_cum{iN}(iT_trial_analysis(3),RT_model(iT_trial_analysis(3))),'k>','markerfacecolor',cmap_c(1,:),'markersize',7);
else
    plot(RT_model(iT_trial_analysis(3)), RespMat_logPR_cum{iN}(iT_trial_analysis(3),RT_model(iT_trial_analysis(3))),'k<','markerfacecolor',cmap_c(3,:),'markersize',7);
end
ylim([-16 16]); xlim([0 POparam.sampleSet(iN)+1]);
line([0 POparam.sampleSet(iN)+1],[0 0],'linestyle','--','color','k'); 
view(90,90); 
set(SP, 'box', 'off', 'TickDir', 'out','FontSize', Common_fontsize,'Xtick',[],'Ytick',[-10 0 10],'xcolor','w');


% 
FigHandle = gcf; 
t = clock;
save_figureToPDF(FigHandle, ['Sampling_algorithm' date() ]);













































    

%% Figure
set(figure(11),'position',[518 443 453 253]); clf; 
set(gcf,'color','w'); 
Common_fontsize = 10; 
SP = subplot('position', [0.1300 0.55 0.7750 0.3]); cla; hold on;
image(1:100,linspace(1,-1,101),flipud(255*(1-result.matPOST'./max(max(result.matPOST))))); hold on;
colormap(gray(256));
line([0 101],[0 0],'linestyle','--','color','k'); 
plot(1:100,result.maxPi*2-1,'-', 'Color',[0 0 0]+.4,'LineWidth',2);

for iT = 1:length(RespMat_logPR_cum{iN}(:,1))
    if RespMat{iN}(iT,2) == 1
        plot(iT, RespMat{iN}(iT,2), 'wo', 'MarkerFaceColor', cmap_c(1,:), 'MarkerSize', 4.5); 
    else
        plot(iT, RespMat{iN}(iT,2), 'wo', 'MarkerFaceColor', cmap_c(3,:), 'MarkerSize', 4.5); 
    end
end

yTickName{1}='ccw';yTickName{2}='cw';
set(gca,'yDir','normal');
ylabel('Percept'); 
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [], 'YTick', [-1 1],'YTickLabel',yTickName, 'FontSize', Common_fontsize,'xcolor','none','ycolor','none')
xlim([0  100]); ylim([-1.3 1.3]); 

xtemp = 73;
iT_trial_analysis = [xtemp xtemp+1 xtemp+2]; 
hold on; 
for iT = 1:length(iT_trial_analysis)
    plot(iT_trial_analysis(iT), 1.3,'kv','color',[0 0 0]+ iT*0.2,'markerfacecolor',[0 0 0]+ 0.2, 'MarkerSize', 4); 
end
% text(iT_trial_analysis(1), 1.56, {'Example trial #1~3'}, 'FontSize', Common_fontsize-2);



% LogPR evidence plot 
SP = subplot('position', [0.1300 0.3 0.7750 0.150]); cla; hold on;

yyaxis right
hold on;
plot(find(vectStable==1), RT_model_plot(vectStable==1),'ko','markerfacecolor','k','markersize',3,'markeredgecolor',[0 0 0]+0.5);
plot(find(vectStable~=1), RT_model_plot(vectStable~=1),'ko','markerfacecolor','k','markersize',3,'markeredgecolor',[0 0 0]+0.5);
set(SP, 'YTick', [],'ycolor','k', 'FontSize', Common_fontsize,'ycolor','none');
ylim([-1.5 24]); xlim([0  100]);
% xlabel('Time (trial)'); 

yyaxis left
% V_x = [0 100 100 0]; 
% V_y = [0 0 24 24]; 
% V = [V_x ; V_y]';
% f = 1:length(V); 
% patch('Faces',f,'Vertices',V,'FaceColor',[0 0 0]+0.9,'edgecolor','none','FaceAlpha',0.8);

for iCase = 1:length(matRun(:,2))
    V_x = [matRun(iCase,4)-0.5 matRun(iCase,4):matRun(iCase,3) matRun(iCase,3)+0.5]; 
    V_y = [0 abs(RespMat_logPR_cum{iN}(matRun(iCase,4):matRun(iCase,3),end))' 0]; 
    V = [V_x ; V_y]'; 
    f = 1:length(V); 
    if matRun(iCase,1) == 1
        patch('Faces',f,'Vertices',V,'FaceColor',cmap_c(1,:),'edgecolor','none','FaceAlpha',0.8);
    else
        patch('Faces',f,'Vertices',V,'FaceColor',cmap_c(3,:),'edgecolor','none','FaceAlpha',0.8);
    end
end
xlim([0  100]); ylim([-1.5 25]);
ylabel({'|logPR| ','& RT '});
set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [0:20:100],'XTicklabel',{'','','','',''}, 'XAxisLocation', 'bottom', 'YTick', [0 20], 'FontSize', Common_fontsize,'ycolor','none');
% set(SP, 'box', 'off', 'TickDir', 'out', 'XTick', [0:20:100], 'XAxisLocation', 'bottom', 'YTick', [0 20], 'FontSize', Common_fontsize,'ycolor','none');



FigHandle = gcf; 
save_figureToPDF(FigHandle, 'Fig2bc');





















% 















































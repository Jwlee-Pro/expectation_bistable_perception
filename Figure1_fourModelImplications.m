%% Figure drawing 

% Figure 1 
% Various implications of expectation
clear all; close all; clc; 

addpath('/Volumes/Data_CSNL/people/JWL/Projects(workspace)/Structure-from-motion/Paper Writing/LeeEtal_ManuSfM/2017_M10_figure/MVPA analysis/ForNN_2018May/Brain2Brain'); 
Setting_param_bhv; 

PE = 0:0.01:1; 
UU.cw = 1-PE; 
UU.ccw = PE; 
EU.cw = -PE.*log2(PE) - (1-PE).*log2(1-PE); 
EU.ccw = -PE.*log2(PE) - (1-PE).*log2(1-PE); 
%% Figure drawing

set(figure(2),'position',[2 654 155 151]); clf; SP = subplot(1,1,1); hold on; 
plot(PE+0.02, PE,'color',cmap_c(3,:),'linewidth',2.8);  
plot(PE, PE,'color',cmap_c(1,:),'linewidth',2.8);
line([0 1],[0 0]+0.5,'linestyle','--','color','k'); 
line([0 0]+0.5,[0 1],'linestyle','--','color','k'); 
xlim([-0.07 1]); ylim([-0.5-0.08 0.5]+0.5); 
% xlabel('Expectation'); 
% ylabel('Expectation signal'); 
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [0 0.5 1], 'XTick', [0 0.5 1], 'YTicklabel', [], 'XTicklabel', [], 'FontSize', 10,'color','none')


% The mean of posterior is changed by expectation 
clear SPost SPost_replay
xaxis = -10:0.001:11; 
for iPE = 1:length(PE)
    % SfM
    Posterior = normpdf(xaxis, PE(iPE)-0.5, 0.4)/sum(normpdf(xaxis, PE(iPE)-0.5, 0.4)); 
    zeroInd = find(xaxis == 0); 
    SPost.cw(iPE) = sum(xaxis(1:(zeroInd-1)).*Posterior(1:(zeroInd-1))); 
    SPost.ccw(iPE) = sum(xaxis((zeroInd+1):length(Posterior)).*Posterior((zeroInd+1):length(Posterior))); 
    
    % Replay
    Posterior = normpdf(xaxis, PE(iPE)-1.0, 0.4)/sum(normpdf(xaxis, PE(iPE)-1.0, 0.4)); 
    zeroInd = find(xaxis == 0); 
    SPost_replay.cw(iPE) = sum(xaxis(1:(zeroInd-1)).*Posterior(1:(zeroInd-1))); 
    Posterior = normpdf(xaxis, PE(iPE), 0.4)/sum(normpdf(xaxis, PE(iPE)-1.0, 0.4)); 
    zeroInd = find(xaxis == 0); 
    SPost_replay.ccw(iPE) = sum(xaxis((zeroInd+1):length(Posterior)).*Posterior((zeroInd+1):length(Posterior))); 
end

set(figure(3),'position',[1 429 155 151]); clf; SP = subplot(1,1,1); hold on; 
plot(1-PE, -SPost.cw+0.5,'color',cmap_c(3,:),'linewidth',2.8);  
plot(1-PE, -SPost.ccw+0.5,'color',cmap_c(1,:),'linewidth',2.8);
line([0 1],[0 0]+0.5,'linestyle','--','color','k'); 
line([0 0]+0.5,[0 1],'linestyle','--','color','k'); 
xlim([-0.07 1]); ylim([-0.5-0.08 0.5]+0.5); 
% xlabel('Expectation'); 
% ylabel('Sampled posterior');
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [0 0.5 1], 'XTick', [0 0.5 1], 'YTicklabel', [], 'XTicklabel', [], 'FontSize', 10,'color','none')

set(figure(30),'position',[1 429 155 151]); clf; SP = subplot(1,1,1); hold on; 
plot(1-PE, -SPost_replay.cw,'color',cmap_c(3,:),'linewidth',2.8);  
plot(1-PE, -SPost_replay.ccw,'color',cmap_c(1,:),'linewidth',2.8);
plot(1-PE, -SPost.cw,'k-.','color',[219 182 104]/255,'linewidth',2.8);  
plot(1-PE, -SPost.ccw,'k-.','color',[111 176 217]/255,'linewidth',2.8);
line([0 1],[0 0],'linestyle','--','color','k'); 
line([0 0]+0.5,[-1.5 1.5],'linestyle','--','color','k'); 
xlim([-0.07 1]); ylim([-1.5 1.5]); 
% xlabel('Expectation'); 
% ylabel('Sampled posterior');
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [0 0.5 1], 'XTick', [0 0.5 1], 'YTicklabel', [], 'XTicklabel', [], 'FontSize', 10,'color','none')


set(figure(5),'position',[2 654 155 151]); clf; SP = subplot(1,1,1); hold on; 
plot(PE+0.02, EU.cw,'color',cmap_c(3,:),'linewidth',2.8);  
plot(PE, EU.ccw,'color',cmap_c(1,:),'linewidth',2.8);
% line([0 1],[0 0]+0.5,'linestyle','--','color','k'); 
line([0 0]+0.5,[0 1],'linestyle','--','color','k'); 
xlim([-0.07 1]); ylim([-0.5-0.08 0.5]+0.5); 
% xlabel('Expectation'); 
% ylabel('Expected uncertainty');
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [0 1], 'XTick', [0 0.5 1], 'YTicklabel', [], 'XTicklabel', [], 'FontSize', 10,'color','none')


set(figure(4),'position',[2 654 155 151]); clf; SP = subplot(1,1,1); hold on; 
plot(PE+0.02, UU.cw,'color',cmap_c(3,:),'linewidth',2.8);  
plot(PE, UU.ccw,'color',cmap_c(1,:),'linewidth',2.8);
% line([0 1],[0 0]+0.5,'linestyle','--','color','k'); 
line([0 0]+0.5,[0 1],'linestyle','--','color','k'); 
xlim([-0.07 1]); ylim([-0.5-0.08 0.5]+0.5); 
% xlabel('Expectation'); 
% ylabel('Unexpected uncertainty');
set(SP, 'box', 'off', 'TickDir','out', 'YTick', [0 1], 'XTick', [0 0.5 1], 'YTicklabel', [], 'XTicklabel', [], 'FontSize', 10,'color','none')





































































































































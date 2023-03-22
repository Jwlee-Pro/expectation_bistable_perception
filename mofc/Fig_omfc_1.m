Percept_all{iSub} = [Percept_all{iSub} matSFM{iS,iSub}.result_rs.Percept]; 
DV_all{iSub} = [DV_all{iSub} matSFM{iS,iSub}.result_rs.maxPi]; 
PE_all{iSub} = [PE_all{iSub} matSFM{iS,iSub}.result_rs.matPrior];         
H_all{iSub} = [H_all{iSub} matSFM{iS,iSub}.result_rs.EU];         
U_all


iSub = 9 ; iS = 1; 
load('mDataMat_all'); 

figure(1); clf; 
SP = subplot(1,1,1); cla; hold on; 
yyaxis left
plot(find(mDataMat.Decision{iSub,iS}==1), ones(1,sum(mDataMat.Decision{iSub,iS}==1)),'wo','markerfacecolor',cmap_c(1,:)); 
plot(find(mDataMat.Decision{iSub,iS}~=1), -ones(1,sum(mDataMat.Decision{iSub,iS}~=1)),'wo','markerfacecolor',cmap_c(3,:)); 
plot(mDataMat.DV{iSub,iS},'m.-'); 
ylabel('DV,decision'); 
set(gca, 'ycolor','k'); 

yyaxis right
plot(mDataMat.PE{iSub,iS},'b.-'); 
ylim([0 1]); 
plot([0 136],[0.5 0.5],'k--'); 
ylabel('E(\omega)'); 
set(gca, 'ycolor','b'); 
xlabel('Time (trials)'); 



prevChoice = mDataMat.Decision{iSub,iS}(1:(end-1)); 
prevPE = mDataMat.PE{iSub,iS}(1:(end-1)); 
currPE = mDataMat.PE{iSub,iS}(2:end);
deltaPE = currPE - prevPE; 
deltaPE = mDataMat.DV{iSub,iS}(1:(end-1)) - mDataMat.DV{iSub,iS}(2:end); 

set(figure(2),'position',[124 471 223 231]); clf; 
SP = subplot(1,1,1); cla; hold on; 
plot(prevPE(prevChoice==1), deltaPE(prevChoice==1),'ro','color',cmap_c(1,:)); 
plot(prevPE(prevChoice~=1), deltaPE(prevChoice~=1),'bo','color',cmap_c(3,:)); 
plot([0 1],[0 0],'k--'); 
xlabel('E_t_-_1(w)'); xlim([0 1]); 
xticks(linspace(0,1,3));
ylabel('\Delta = E_t(w) - E_t_-_1(w)'); 








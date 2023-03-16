function vectT = PrefTagging_from_sfm_FFFF(TS, paramSFM, iSub)
% TS: replay TS [iVox x iT] 

TR = 2.4; 
HemodynamicDelay = 2*TR; 

truth_fMRI_all = [] ; 
% Pre-process behavioral data 
	for iS =1:4
		strInd = find(paramSFM.track{iS}(:,2) == 0); 
		if ~isempty(strInd)
			for iSt = 1:length(strInd)
				paramSFM.track{iS}(strInd(iSt),2) = paramSFM.track{iS}(strInd(iSt)-1,2) + paramSFM.tBlank + paramSFM.tDisplay; 
			end
		end
			

		truth = paramSFM.track{iS}(:,1); 
		truth_time = paramSFM.track{iS}(:,2); 


    

		% Assigning perception to each TR
		time_fMRI = (TR:TR:TR*(length(TS{iS}(1,:))));
		minT = min(truth_time + HemodynamicDelay); maxT = max(truth_time + HemodynamicDelay); 
		minT_fMRI = find(time_fMRI > minT,1,'first'); maxT_fMRI = find(time_fMRI < maxT, 1, 'last');
		
		% Undersample this (since trial length < sampling period) 
		% A little trick for SFM because of "ambiguous" condition 
		truth_fMRI = zeros(1,length(TS{1}(1,:))-(HemodynamicDelay/TR)); 
		for iT = 1:length(TS{1}(1,:))-(HemodynamicDelay/TR)
			timeInd = time_fMRI(iT); % after correcting for hemodynamic delay 
			pInd_bf = find(timeInd > truth_time,1,'last'); 
			if (truth(pInd_bf) == 1) & (truth(pInd_bf+1) == 1)
				truth_fMRI(iT) = 1; % left
			elseif (truth(pInd_bf) == -1) & (truth(pInd_bf+1) == -1)
				truth_fMRI(iT) = 2; % right
			elseif (truth(pInd_bf) == 1) & (truth(pInd_bf+1) == -1)
				truth_fMRI(iT) = 1; % follow the pInd_bf. 
			elseif (truth(pInd_bf) == -1) & (truth(pInd_bf+1) == 1)
				truth_fMRI(iT) = 2; % follow the pInd_bf. 
			end
		end
		
		truth_fMRI_set{iS} = truth_fMRI;  
        truth_fMRI_all = [truth_fMRI_all truth_fMRI_set{iS}];
    end
	
    L_ind = find(truth_fMRI_all == 1); % L perception
    R_ind = find(truth_fMRI_all == 2); % R perception

		
    % Generate the TS 
    TS_all = []; 
    for iS = 1:4
        TS_all = [TS_all TS{iS}(:, 1 + (HemodynamicDelay/TR):end)] ; 
    end
    
    % Response space (RV)
    x = -5:0.01:5;

    
    % Decode information from population activity 
    for iVox = 1:length(TS_all(:,1)) 
        rL = TS_all(iVox, L_ind);
        rR = TS_all(iVox, R_ind);
        
        f_L = exp(-(x-mean(rL)).^2./(2*std(rL)^2))./(std(rL)*sqrt(2*pi)); % Pr(r|p='L')
        f_R = exp(-(x-mean(rR)).^2./(2*std(rR)^2))./(std(rR)*sqrt(2*pi)); % Pr(r|p='R')
    
        
        % Pr(p='L'|r)
        um_L = f_L * length(rL)/(length(rL)+length(rR)); 
        P_r = (f_L * length(rL)/(length(rL)+length(rR))) + (f_R * length(rR)/(length(rL)+length(rR))); 
        f_L_posterior(iVox,:) = um_L./P_r;


        nm = ((length(L_ind)-1)*(std(rL).^2) + (length(R_ind)-1)*(std(rR).^2));
        dm = (length(L_ind) + length(R_ind) - 2);
        se_pooled = sqrt(nm./dm);
        ft = 1/(sqrt((1/length(L_ind)) + (1/length(R_ind))));

        diff_mean = mean(rL)-mean(rR); 

%         vectT(iVox) = ft*(diff_mean./se_pooled);
        vectT(iVox) = diff_mean;
        
        rAll = [rL rR] ; 
    end

    


    
% Find voxels with max preference (for figure) 
    maxVoxInd = find(vectT == max(vectT)); 
    minVoxInd = find(vectT == min(vectT)); 
    zeroVoxInd = find(abs(vectT) == min(abs(vectT))); 



    
% Figure 
% subplot(1); bias
% subplot(2); tagged histogram
% subplot(3); preference tagging result

% figure(201); 
% subplot(2,14,iSub); 
% 
% subplot(2,14,14 + iSub); 
% hold on; 
% hist(vectT); 
% set(get(gca,'child'),'FaceColor',[0 0 0]+0.7,'EdgeColor','w');

% subplot(4,14,2*14 + iSub); 
% hist(vectT); 
% set(get(gca,'child'),'FaceColor',[0 0 0]+0.7,'EdgeColor','w');




% figure();  
% clf;
% 
% subplot(1,3,1); % right like
% hist(TS(maxVoxInd,L_ind)); hold on; 
% hist(TS(maxVoxInd,R_ind));
% h = findobj(gca,'Type','patch');
% display(h);
% % set(h(1),'FaceColor','none','EdgeColor','r','facealpha',0.5); 
% % set(h(2),'FaceColor','none','EdgeColor','b','facealpha',0.5);
% set(h(1),'FaceColor','b','EdgeColor','w','facealpha',0.3); 
% set(h(2),'FaceColor','r','EdgeColor','w','facealpha',0.3); 
% title('Left-preferred','fontweight','bold','fontsize',12);
% xlabel('BOLD amplitude (%)'); ylabel('# of samples');    
% legend('Activity to Left p','Activity to Right p'); 
% line([0 0]+mean(TS(maxVoxInd,L_ind)),[0 18],'color','r','linewidth',1.5); 
% line([0 0]+mean(TS(maxVoxInd,R_ind)),[0 18],'color','b','linewidth',1.5); 
% xlim([-3 3]); ylim([0 18]); 
% 
% subplot(1,3,2); % same same
% hist(TS(zeroVoxInd,L_ind)); hold on; 
% hist(TS(zeroVoxInd,R_ind));
% h = findobj(gca,'Type','patch');
% display(h);
% % set(h(1),'FaceColor','none','EdgeColor','r','facealpha',0.5); 
% % set(h(2),'FaceColor','none','EdgeColor','b','facealpha',0.5);
% set(h(1),'FaceColor','b','EdgeColor','w','facealpha',0.3); 
% set(h(2),'FaceColor','r','EdgeColor','w','facealpha',0.3); 
% title('No preference','fontweight','bold','fontsize',12);
% xlabel('BOLD amplitude (%)'); ylabel('# of samples');    
% xlim([-3 3]); ylim([0 18]); 
% line([0 0]+mean(TS(zeroVoxInd,L_ind)),[0 18],'color','r','linewidth',1.5); 
% line([0 0]+mean(TS(zeroVoxInd,R_ind)),[0 18],'color','b','linewidth',1.5); 
% 
% subplot(1,3,3); % left like
% hist(TS(minVoxInd,L_ind)); hold on; 
% hist(TS(minVoxInd,R_ind));
% h = findobj(gca,'Type','patch');
% display(h);
% % set(h(1),'FaceColor','none','EdgeColor','r','facealpha',0.5); 
% % set(h(2),'FaceColor','none','EdgeColor','b','facealpha',0.5);
% set(h(1),'FaceColor','b','EdgeColor','w','facealpha',0.3); 
% set(h(2),'FaceColor','r','EdgeColor','w','facealpha',0.3); 
% title('Right-preferred','fontweight','bold','fontsize',12);
% xlabel('BOLD amplitude (%)'); ylabel('# of samples');    
% xlim([-3 3]); ylim([0 18]); 
% line([0 0]+mean(TS(minVoxInd,L_ind)),[0 18],'color','r','linewidth',1.5); 
% line([0 0]+mean(TS(minVoxInd,R_ind)),[0 18],'color','b','linewidth',1.5); 


end

function [xTime, DM, DM2] = GenSurvivalRate_merging(track, tTrial, tCrit)
xInd = 1:floor(tCrit/tTrial); 
xTime = xInd*tTrial; 

DM  = nan(length(xInd), 2); 
DM2 = nan(length(xInd), 2); 

for wInd = 1:length(xInd)
    ccnt_s   = 0; 
    ccnt     = 0; 
    ccnt_s2  = 0;
    ccnt2    = 0;
    
    if ~isempty(track)
        for iInd = 1:(length(track(:,1))-wInd)
            if track(iInd,1) == track(iInd + wInd,1)
                ccnt_s = ccnt_s + (track(iInd,1) == 1); 
                ccnt_s2 = ccnt_s2 + (track(iInd,1) == -1); 
            end
            ccnt = ccnt + (track(iInd,1) == 1);
            ccnt2 = ccnt2 + (track(iInd,1) == -1);
        end
    end

    DM(wInd, 1) = ccnt_s; 
    DM(wInd, 2) = ccnt; 
    DM2(wInd, 1) = ccnt_s2; 
    DM2(wInd, 2) = ccnt2; 
end

% % Calculate survival rate
% SR_temp_cw   = DM(:, 1)./DM(:, 2); 
% SR_temp_ccw  = DM2(:, 1)./DM2(:, 2); 
% 
% DataAll = DM + DM2;
% 
% SR{1} = SR_temp_cw;
% SR{2} = SR_temp_ccw;
% % Total perception
% SR{3} = DataAll(:, 1)./DataAll(:, 2); 

    

    
    function [xInd, SR] = GenSurvivalRate_trial(track, tTrialCrit)
xInd = 1:tTrialCrit; 

for wInd = 1:length(xInd)
    ccnt_s   = 0; 
    ccnt     = 0; 
    ccnt_s2  = 0;
    ccnt2    = 0;
    for iS = 1:length(track)
        for iInd = 1:length(track{iS}(:,1))-wInd
            if track{iS}(iInd,1) == 1 % cw
                if track{iS}(iInd,1) == track{iS}(iInd + wInd,1)
                    ccnt_s = ccnt_s + 1; 
                end
                ccnt = ccnt + 1; 
            elseif track{iS}(iInd,1) == -1 % ccw
                if track{iS}(iInd,1) == track{iS}(iInd + wInd,1)
                    ccnt_s2 = ccnt_s2 + 1; 
                end
                ccnt2 = ccnt2 + 1; 
            end
        end
    end
    DataMat(wInd, 1) = ccnt_s; 
    DataMat(wInd, 2) = ccnt; 
    DataMat2(wInd, 1) = ccnt_s2; 
    DataMat2(wInd, 2) = ccnt2; 
end

% Calculate survival rate
SR_temp_cw   = DataMat(:, 1)./DataMat(:, 2); 
SR_temp_ccw  = DataMat2(:, 1)./DataMat2(:, 2); 

DataAll = DataMat + DataMat2;

SR{1} = SR_temp_cw;
SR{2} = SR_temp_ccw;
% Total perception
SR{3} = DataAll(:, 1)./DataAll(:, 2); 

    
    
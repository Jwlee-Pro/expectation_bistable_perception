function [xTime, TR] = GenCumTransitionRate_sess(track, tTrial, tCrit)
xInd = 1:floor(tCrit/tTrial); 
xTime = xInd*tTrial; 

for wInd = 1:length(xInd)
    ccnt_s   = 0; 
    ccnt     = 0; 
    ccnt_s2  = 0;
    ccnt2    = 0;
    for iS = 1:length(track)
        if ~isempty(track)
            DiffMat = [0; abs(diff(track{iS}(:,1)))]; 
            for iInd = 1:(length(track{iS}(:,1))-wInd)
                if sum(DiffMat(iInd:(iInd + wInd))) ~= 0 
                    ccnt_s  = ccnt_s  + (track{iS}(iInd,1) == 1);  % cw
                    ccnt_s2 = ccnt_s2 + (track{iS}(iInd,1) == -1); % ccw
                end
                ccnt = ccnt + (track{iS}(iInd,1) == 1);
                ccnt2 = ccnt2 + (track{iS}(iInd,1) == -1);
            end
        end
    end
    TM(wInd, 1) = ccnt_s; 
    TM(wInd, 2) = ccnt; 
    TM2(wInd, 1) = ccnt_s2; 
    TM2(wInd, 2) = ccnt2; 
end
            

% Calculate transition rate
SR_temp_cw   = TM(:, 1)./TM(:, 2); 
SR_temp_ccw  = TM2(:, 1)./TM2(:, 2); 

DataAll = TM + TM2;

TR{1} = SR_temp_cw;
TR{2} = SR_temp_ccw;
% Total perception
TR{3} = DataAll(:, 1)./DataAll(:, 2); 

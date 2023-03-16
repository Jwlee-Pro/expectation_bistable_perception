
%% Behavioral data preprocessing and plotting

function [paramOutput] = behPreprocess_fmri_old(paramSFM)
    for iS = 1:paramSFM.nSession
        strInd = find(paramSFM.track{1,iS}(:,2) == 0); 
        if ~isempty(strInd)
            for iSt = 1:length(strInd)
                if strInd(iSt) ~= 1
                    paramSFM.track{1,iS}(strInd(iSt),2) = paramSFM.track{1,iS}(strInd(iSt)-1,2) + paramSFM.tBlank + paramSFM.tDisplay; 
                else
                    paramSFM.track{1,iS}(strInd(iSt),2) = 1; 
                end
            end
        end
    end
    paramOutput = paramSFM; 
end
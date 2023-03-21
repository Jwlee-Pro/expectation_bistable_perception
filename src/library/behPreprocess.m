
%% Behavioral data preprocessing and plotting

function [paramOutput] = behPreprocess(paramSFM)
    strInd = find(paramSFM.track{1,1}(:,2) == 0); 
    if ~isempty(strInd)
        for iSt = 1:length(strInd)
            if strInd(iSt) ~= 1
                paramSFM.track{1,1}(strInd(iSt),2) = paramSFM.track{1,1}(strInd(iSt)-1,2) + paramSFM.tBlank + paramSFM.tDisplay; 
            else
                paramSFM.track{1,1}(strInd(iSt),2) = 1; 
            end
        end
    end
    paramOutput = paramSFM; 
end
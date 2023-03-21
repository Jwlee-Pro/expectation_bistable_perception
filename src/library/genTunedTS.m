function TSout = genTunedTS(TS, RegTS, PrefInd)
% Input: TS [nVox x nT], PrefInd [1 x nVox]
% Output: TSout [1 x nT]
    TSresi = [] ; 
    % Regressing out by mean 
    for iVox = 1:length(TS(:,1))
%         matTS = [ones(length(TS(1,:)),1) RegTS' (RegTS.^2)']; 
        matTS = [ones(length(TS(1,:)),1) RegTS']; 
        [b, bint, resi, rint, rstat] = regress(TS(iVox,:)', matTS);
%         resi = TS(iVox,:)-RegTS;
        TSresi(iVox,:) = resi; 
    end
    
    % Calculate slope 
    for iT = 1:length(TS(1,:))
        temp = polyfit(PrefInd, (TSresi(:,iT)'), 1);

%          temp = polyfit(PrefInd, (TSresi(:,iT)'), 1);
%          temp = polyfit(PrefInd, abs(PrefInd).*TSresi(:,iT)', 1);
%          temp = polyfit(PrefInd, 2*(TSresi(:,iT)'>0) -1, 1);

%            temp = mean(TSresi(find(PrefInd>0),iT)) - mean(TSresi(find(PrefInd<0),iT)) ; 
                      


%            if (sum(PrefInd>0) > 0)
%                val1 = mean(TSresi(find(PrefInd>0),iT));
%            else
%                val1 = 0; 
%            end
%            if (sum(PrefInd<0) > 0)
%                val2 = mean(TSresi(find(PrefInd<0),iT));
%            else
%                val2 = 0; 
%            end
%            temp = val1 - val2 ; 
           

           
           
%            absPref = abs(PrefInd); 
%            if (sum(PrefInd>0) > 0)
%                val1 = mean(absPref(find(PrefInd>0))'.*TSresi(find(PrefInd>0),iT));
%            else
%                val1 = 0; 
%            end
%            if (sum(PrefInd<0) > 0)
%                val2 = mean(absPref(find(PrefInd<0))'.*TSresi(find(PrefInd<0),iT));
%            else
%                val2 = 0; 
%            end
%            temp = val1 - val2 ; 


%            temp = abs(mean(TSresi(find(PrefInd>0),iT)))/abs(mean(TSresi(find(PrefInd<0),iT))) ; 

         TSout(1,iT) = temp(1);     
    end
end


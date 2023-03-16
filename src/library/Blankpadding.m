function pupil_padded = Blankpadding(pupildata)
    Ind_bfnan = find(diff(isnan(pupildata))==1);
    Ind_afnan = find(diff(isnan(pupildata))==-1)+1;
    if isnan(pupildata(end)) == 1
        Ind_afnan = [Ind_afnan; length(pupildata)];
    end
    nNan = (min(length(Ind_bfnan), length(Ind_afnan)));
    if length(Ind_bfnan) > length(Ind_afnan) 
        for iBlink = 1:nNan
            pupildata(Ind_bfnan(iBlink)+1:Ind_afnan(iBlink)) = linspace(pupildata(Ind_bfnan(iBlink)), pupildata(Ind_afnan(iBlink)), Ind_afnan(iBlink)-Ind_bfnan(iBlink)); 
        end
        pupildata(Ind_bfnan(iBlink+1):end) = pupildata(Ind_bfnan(iBlink+1)); 
    elseif length(Ind_bfnan) == length(Ind_afnan)
        for iBlink = 1:nNan
            pupildata(Ind_bfnan(iBlink)+1:Ind_afnan(iBlink)) = linspace(pupildata(Ind_bfnan(iBlink)), pupildata(Ind_afnan(iBlink)), Ind_afnan(iBlink)-Ind_bfnan(iBlink));
        end
    elseif length(Ind_bfnan) < length(Ind_afnan)
        for iBlink = 1:nNan
            pupildata(Ind_bfnan(iBlink)+1:Ind_afnan(iBlink+1)) = linspace(pupildata(Ind_bfnan(iBlink)), pupildata(Ind_afnan(iBlink+1)), Ind_afnan(iBlink+1)-Ind_bfnan(iBlink)); 
        end
        pupildata(1:Ind_afnan(1)) = pupildata(Ind_afnan(1)+1); 
    end
    if isnan(pupildata(end))
        pupildata(find(isnan(pupildata),1,'first'):length(pupildata)) = pupildata(find(isnan(pupildata),1,'first')-1); 
    end
    pupil_padded = pupildata; 
end
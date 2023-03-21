function TS_noAmb = Amb2twoAFC(TS) 
    % input: TS [nT x 2] 
    % output: TS_noAmb [nT x 2]

    TS_noAmb = TS; 

    AmbInd = find(TS(:,1) == 3); 


    for iA = 1:length(AmbInd)
        if (AmbInd(iA) == 1) & (TS(AmbInd(iA)+1,1) ~= 3)
            TS_noAmb(AmbInd(iA),1) = -TS_noAmb(AmbInd(iA)+1,1); 
        elseif (AmbInd(iA) == 1) & (TS_noAmb(AmbInd(iA)+1,1) == 3)
            TS_noAmb(AmbInd(iA),1) = 1; % just assign 'Left'
        else
            TS_noAmb(AmbInd(iA),1) = -TS_noAmb(AmbInd(iA)-1,1); 
        end
    end

    
   NoRespInd = find(TS(:,1) == 0); 
   for iN = 1:length(NoRespInd)
       if NoRespInd(iN) == 1
           TS_noAmb(NoRespInd(iN),1) = 1; 
       else
           TS_noAmb(NoRespInd(iN),1) = TS_noAmb(NoRespInd(iN)-1,1); 
       end
   end
   

end





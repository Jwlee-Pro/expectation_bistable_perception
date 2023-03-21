function temp = plot_StackHist(DataSet1, DataSet2)

    DataHist = [DataSet1 DataSet2]; 
    minHist = min(DataHist); maxHist = max(DataHist); 
    BinRange = linspace(minHist, maxHist, 50); 
    % center 
    for iB = 1:length(BinRange)-1
        BinCenter(iB) = (BinRange(iB)+BinRange(iB+1))/2; 
    end
    BinCenter_sided = [(2*BinCenter(1) - (BinCenter(2))) BinCenter (2*BinCenter(end) - (BinCenter(end-1)))];

    figure(100); 

    hInfo = histogram(DataSet1,BinRange); 
    temp.BinVal_1 = hInfo.Values; 
    hInfo = histogram(DataSet2,BinRange); 
    temp.BinVal_2 = hInfo.Values; 
    temp.BinCenter = BinCenter; 
    temp.minHist = minHist; 
    temp.maxHist = maxHist; 
    

end

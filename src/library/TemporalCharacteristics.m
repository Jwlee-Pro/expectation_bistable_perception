function [Peakloc, WidthLong, WidthRange] = TemporalCharacteristics(TSx, TS_time)
global gParam time_fmri session

BaseInd = find(abs(TS_time+(gParam.TR*gParam.HD)) == min(abs(TS_time+(gParam.TR*gParam.HD)))); 
MaxInd = find(TSx == max(TSx)); 
peakAmp = (max(TSx)-TSx(BaseInd));
TS_time((abs(TSx-(2*peakAmp/3 + TSx(BaseInd))) < 0.1));

% after point
tempx = find(abs(TSx-(2*peakAmp/3 + TSx(BaseInd))) < 0.1);
tempx1 = tempx(tempx > MaxInd); 
finalInd = tempx1(abs(TSx(tempx1) -(2*peakAmp/3 + TSx(BaseInd))) == min(abs(TSx(tempx1) -(2*peakAmp/3 + TSx(BaseInd))))); 

% before point 
tempx1 = tempx((tempx < MaxInd) & (tempx > BaseInd)); 
preInd = tempx1(abs(TSx(tempx1) -(2*peakAmp/3 + TSx(BaseInd))) == min(abs(TSx(tempx1) -(2*peakAmp/3 + TSx(BaseInd))))); 

Peakloc = TS_time(MaxInd);
WidthLong = TS_time(finalInd) - TS_time(preInd);

WidthRange = [TS_time(preInd) TS_time(finalInd)];

function ctSeries = removeLowFqDrift(tSeries, period)
% function tSeries = removeBaseline(tSeries, period);
%
% Use multiple boxcar smoothing operations to remove low-frequency baseline
% drift from the input time series (tSeries) using the input period in
% FRAMES (not seconds!). The input is assumed to be  two-dimensional, with
% time as the low-order dimension.
%
% DBR  5/00

% djh & dbr, 2/2001. 
% Fixed numIterations=2, no longer a variable input argument.
%
[nFrame, nEl] = size(tSeries);
%
% Create boxcar kernel:
kernel = ones([period 1]) / period;
numIterations = 2;

% Initialize the baseline array to the time with 1-period
% padding at beginning and end:
ntPoints = size(tSeries,1);
nSeries = size(tSeries,2);
mValues = mean(tSeries);
nBLine = ntPoints + 2*period;
bLine = zeros(nBLine,nSeries);
firstTrialMean = mean(tSeries(1:period,:));

for frame = 1:period
    bLine(frame,:) = firstTrialMean; 
end
bLine(period+1:period+ntPoints,:) = tSeries;
lastTrialMean = mean(tSeries(ntPoints-period+1:ntPoints,:));
for frame = period+ntPoints+1:nBLine
    bLine(frame,:) = lastTrialMean; 
end

% Define indices for post-smoothing array "trim":
addPts = numIterations * (period - 1);
start = floor(addPts/2) + 1;
stop = nBLine + floor(addPts/2);

% Smoothing loop -- convolve with boxcar, then "trim" array:
for i=1:numIterations
    bLine = conv2(bLine, kernel); 
end
bLine = bLine(start:stop, :);

% Remove baseline from time series:
lowFrComp = bLine(period+1:period+ntPoints,:) ;
%
lowFrComp = lowFrComp - ones(nFrame,1) * mean(tSeries,1) ;
%
lowFrComp = lowFrComp - mean(lowFrComp,1);
%
ctSeries = tSeries - lowFrComp;

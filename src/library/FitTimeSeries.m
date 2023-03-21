function [paramSet, FittedTS] = FitTimeSeries(target, seed, modelNum)
% modelNum1 = surprise
% modelNum2 = entropy
if modelNum < 3
    paramModel1 = [10000, 0, 0]; % surprise model (weight, baseline)
    [paramout_temp, fval_temp] = fminsearchbnd( @(param) fitModel1(param, target, seed, modelNum), paramModel1, [0 -10 -1] ,[Inf 10 1]);
    % Once more to get higher accuracy
    [paramSet, fval_temp] = fminsearchbnd( @(param) fitModel1(param, target, seed, modelNum), paramout_temp, [0 -10 -1] ,[Inf 10 1]);
else
    paramModel1 = [10000, 10000, 0, 0]; % surprise model (weight, baseline)
    [paramout_temp, fval_temp] = fminsearchbnd( @(param) fitModel1(param, target, seed, modelNum), paramModel1, [0 0 -10 -1] ,[Inf Inf 10 1]);
    % Once more to get higher accuracy
    [paramSet, fval_temp] = fminsearchbnd( @(param) fitModel1(param, target, seed, modelNum), paramout_temp, [0 0 -10 -1] ,[Inf Inf 10 1]);
end
% Model structure = polynomial (1st order) with lineare slope 
if modelNum < 3
    FittedTS = paramSet(1)*seed{modelNum} + paramSet(2) +  paramSet(3)*(1:length(seed{modelNum})); 
else
    FittedTS = paramSet(1)*seed{1} + paramSet(2)*seed{2} + paramSet(3) +  paramSet(4)*(1:length(seed{1})); 
end










function [TS_ppc_fin, TS_ppc_th_fin] = LOOcrossValidate(TS_cat, map_covert2TSorder, nGroup)

global gParam time_fmri


% Calculate CP, leave-one-out crossvalidation, 10 groups 
percept_ind_sum = [ones(1,length(TS_cat.cw(:,1))) -ones(1,length(TS_cat.ccw(:,1)))]';
TS_sum = [TS_cat.cw; TS_cat.ccw]; 
map_rand = randperm(length(TS_sum(:,1))); 
[~,map_reverse] = sort(map_rand); 

% Randomized percept and TS 
percept_ind_sum_rand = percept_ind_sum(map_rand); 
TS_sum_rand = TS_sum(map_rand,:); 

nTperGroup = floor(length(TS_sum(:,1))/nGroup); 
nVmat = zeros(1,nGroup) + nTperGroup; nVmat_cum = cumsum(nVmat) - nVmat(1); 
nMod = mod(length(TS_sum(:,1)), nGroup); 
tempG = randperm(nGroup); 
for iM = 1:nMod
    nVmat(tempG(iM)) = nVmat(tempG(iM)) + 1; 
end

TSdec_group = cell(1,5); % there are 5 types of decoder
TSdec_th_group = cell(1,5); 

for iG = 1:nGroup 
    indTest = zeros(1,length(TS_sum(:,1))); 
    if iG ~= nGroup
        indTest((nVmat_cum(iG)+1):nVmat_cum(iG+1)) = 1; 
    else
        indTest((nVmat_cum(iG)+1):end) = 1; 
    end

    matTest.p  = percept_ind_sum_rand(indTest==1); 
    matTest.ts = TS_sum_rand(indTest==1,:); 

    matTrain.p  = percept_ind_sum_rand(indTest~=1); 
    matTrain.ts = TS_sum_rand(indTest~=1,:);

    [tempX, tempY, Perform_group(1, iG)] = genChoiceSignalsimple(matTest, matTrain); 

    
    
    
    TSdec_group{1} = [TSdec_group{1} tempX.meanDiff]; 
    TSdec_group{2} = [TSdec_group{2} tempX.slope]; 
    TSdec_group{3} = [TSdec_group{3} tempX.weightedDiff]; 
    TSdec_group{4} = [TSdec_group{4} tempX.bayesDec]; 
    TSdec_group{5} = [TSdec_group{5} tempX.svm]; 

    TSdec_th_group{1} = [TSdec_th_group{1} tempY.meanDiff]; 
    TSdec_th_group{2} = [TSdec_th_group{2} tempY.slope]; 
    TSdec_th_group{3} = [TSdec_th_group{3} tempY.weightedDiff]; 
    TSdec_th_group{4} = [TSdec_th_group{4} tempY.bayesDec]; 
    TSdec_th_group{5} = [TSdec_th_group{5} tempY.svm]; 
end


for iDecMethod = 1:length(TSdec_th_group)
    % Reorder data
    Tempxx1 = TSdec_group{iDecMethod}(map_reverse); 
    Tempxx2 = TSdec_th_group{iDecMethod}(map_reverse); 
    X1 = Tempxx1(map_covert2TSorder); 
    X2 = Tempxx2(map_covert2TSorder); 
    for iS = 1:length(X1)/(125-gParam.HD)
        TS_ppc_fin{iDecMethod}{1,iS} = X1((((125-gParam.HD)*(iS-1)+1) : (125-gParam.HD)*(iS))); 
        TS_ppc_th_fin{iDecMethod}{1,iS} = X2((((125-gParam.HD)*(iS-1)+1) : (125-gParam.HD)*(iS))); 
	end
end





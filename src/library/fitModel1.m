function errVal = fitModel1(param, target, seed, modelType)
    if modelType == 1
        temp = param(1)*seed{1}(1:end) + param(2) + param(3)*(1:length(seed{1}(1:end))); 
    elseif modelType == 2
        temp = param(1)*seed{2}(1:end) + param(2) + param(3)*(1:length(seed{2}(1:end))); 
    elseif modelType == 3
        temp = param(1)*seed{1}(1:end) + param(2)*seed{2}(1:end) + param(3) + param(4)*(1:length(seed{1}(1:end)));
    end
    errVal = sum((temp - target).^2); 
end
    



function errVal = fitModel_HRFderiv(param, target, seed)
    temp = param(1)*seed{1}(1:end) + param(2)*seed{2}(1:end) + param(3) + param(4)*(1:length(seed{1}(1:end)));
    errVal = sum((temp - target).^2); 
end
    



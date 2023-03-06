function errVal = fitExpo(param, target, seed)
    xx = param(1)*exp(param(2)*seed{1}); 
    errVal = sum((xx - target{1}).^2) + sum((xx - target{2}).^2); 
end
    



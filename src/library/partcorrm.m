function r = partcorrm(x,y,Z,type)
%
% Data should be represented column vectors in x, y, Z.
%

if ~exist('type','var')
    type = 'Pearson';
end

if strcmpi(type,'spearman')
    [~,X]=ismember(X,sort(X,'ascend'));
    [~,Y]=ismember(Y,sort(Y,'ascend'));
%     [~,Z]=ismember(Z,sort(Z,'ascend'));
    
    for iC = 1:size(Z,2)
        [~,Z(:,iC)]=ismember(Z(:,iC),sort(Z(:,iC),'ascend'));
    end
end

% 1) regress out Z from Y
Z_new     = [ones(size(y,1),1), Z]; % make a regressor for Z
[~,~,r_y] = regress(y, Z_new);      % residual from first regression


% 2) regress X with residual
r_new       = [ones(size(r_y,1),1), r_y];   % make a regressor for r_y
[b_x,~,r_x] = regress(x, r_new);


% 3) explained variance (r-squared for the second regression)
r2  = 1 - var(r_x)/var(x);
r   = sign(b_x(2)) * sqrt(r2);


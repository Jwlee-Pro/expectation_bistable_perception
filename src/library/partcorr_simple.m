function [coef,sigout] = partcorr_simple(Y,X,Z,type)
%PARTCORR Part or semipartial correlation coefficient
%   RHO = PARTCORR(Y,X,Z) returns the sample linear part correlation cofficient
%   between Y and X while controlling Z. Y is the response (dependent) variable,
%   X is the explanatory (independent) variable, and Z is the explanatory variable 
%   to be parted out (control).
%
%   Part correlation is the correlation between Y and the residual of X from the linear
%   regression of X by Z.
%
%   Example (from Abdi, H. (2007). Part and partial correlations. In N.J. Salkind (Ed.): Encyclopedia of Measurement and Statistics. Thousand Oaks (CA): Sage. pp. 736-740.)
%         y = [14 23 30 50 39 67]';   % memory span (DV)
%         x = [ 4  4  7  7 10 10]';   % age
%         t = [ 1  2  2  4  3  6]';   % speech rate
%
%         X     = [ones(size(x)), x, t];
%         b     = (X' * X) \ X' * y;
%         y_hat = X * b;
%
%
%         R2.y_xt = 1 - var(y_hat - y) / var(y);    % Partial regression coefficient, squared part correlation (increment in r^2 by adding X in the model)
%
%         r.x_t = corr(x,t);
%         r.y_x = corr(y,x);
%         r.y_t = corr(y,t);
%
%         r2.yx_t = R2.y_xt - r.y_t ^ 2 % Part correlation coefficient. This should be approx. 0.008528
%         r2.yt_x = R2.y_xt - r.y_x ^ 2 % Part correlation coefficient. This should be approx. 0.342072
%
%
%     Then,
%         r2.yx_t should match partcorr(y,x,t)^2
%         r2.yt_x should match partcorr(y,t,x)^2
%
%
%   2014/10/05 RYU Exmaple for validation added, F-stat added
%   2015/06/03 RYU Part rank correlation added
%

if ~exist('type','var')
    type = 'Pearson';
end

if strcmpi(type,'spearman')
    [~,X]=ismember(X,sort(X,'ascend'));
    [~,Y]=ismember(Y,sort(Y,'ascend'));
    [~,Z]=ismember(Z,sort(Z,'ascend'));
end

r.YX = corr(Y,X,'type',type);
r.YZ = corr(Y,Z,'type',type);
r.XZ = corr(X,Z,'type',type);

% replace NaN to zero (e.g. partcorr(ones(10,1),randn(10,1),randn(10,1)) or partialcorr(randn(10,1),randn(10,1),ones(10,1)))
if(isnan(r.YX))
    r.YX = zeros(size(r.YX));
end

if(isnan(r.YZ))
    r.YZ = zeros(size(r.YZ));
end

if(isnan(r.XZ))
    r.XZ = zeros(size(r.XZ));
end

%
coef = (r.YX - r.YZ .* r.XZ) ./ sqrt(1 - r.XZ .^2);
N    = length(Y);
K    = 2;                       % Number of independent variables. for X, Z ==> K=2

stat.R2 = coef ^ 2;             % Variance in Y explained by X, while controlling Z. Increment in explained variance of Y by adding X

if strcmpi(type,'spearman')
    stat.t = coef .* sqrt((N-1)./(1-coef.^2));
else
    stat.F = stat.R2 / (1 - stat.R2 - r.YZ ^ 2) * (N - K - 1);
    stat.t = sqrt(stat.F);
end

stat.p = 1 - tcdf(stat.t, N - K - 1);    % Significance of increment in explained variance by X
sigout = stat.p; 
% Z    = atanh(coef);         % Fisher r-to-z transform
% % Z*sqrt(N-4) should be a normal distribution of stdev 1 with zero mean.
% 
% stat.pval = 1-tcdf(Z*sqrt(N-4),N);	% H0: different from zero
% stat.R2   = coef .^ 2;      % Explained variance of Y by X (vice versa), while controlling Z
% % stat.F    = ;           % Variance reduction analysis
% % stat.F_pval = ;         % significance of variance reduction



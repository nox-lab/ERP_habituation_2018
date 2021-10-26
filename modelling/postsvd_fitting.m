function [fitresult, gof, bic] = postsvd_fitting(V)

%% ERP modelling code relative to:
% Mancini F, Pepe A, Bernacchia, A, Di Stefano G, Mouraux A, Iannetti GD. (2017)
% Characterising the short-term habituation of event-related evoked potentials
% eNeuro

% written in Matlab R2016b by F Mancini, fm456@cam.ac.uk
% This code requires the Curve Fitting Toolbox in Matlab

% Output:
%      fitresult : a structure of fit objects representing the fits.
%      gof : a structure with goodness-of fit info.
%      bic: Bayesian Information Criterion, used for model comparison

x = double([1:size(V,1)]');     % this is the trial number

for rank_ord = 1 : size(V,2)
    
    [fitresult{1,rank_ord}, gof{1,rank_ord}, bic(1,rank_ord)] = eq_one_fit(x, double(V(:,rank_ord)));    % y = a+b/x
    [fitresult{2,rank_ord}, gof{2,rank_ord}, bic(2,rank_ord)] = eq_two_fit(x, double(V(:,rank_ord)));    %  y = a+b/(x^c)
    [fitresult{3,rank_ord}, gof{3,rank_ord}, bic(3,rank_ord)] = eq_three_fit(x, double(V(:,rank_ord)));  %  y = a*exp(-b*x)+c
    [fitresult{4,rank_ord}, bic(4,rank_ord)] = eq_four_fit_ab(x, double(V(:,rank_ord)));                 % y = mean(y) --> no habituation
    
end

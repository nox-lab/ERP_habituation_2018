function [fitresult, bic] = eq_four_fit_ab(x, y)
%% ERP modelling code relative to:
% Mancini F, Pepe A, Bernacchia, A, Di Stefano G, Mouraux A, Iannetti GD. (2017)
% Characterising the short-term habituation of event-related evoked potentials
% bioRxiv doi:10.1101/153387

% written in Matlab R2016b by F Mancini, fm456@cam.ac.uk
% This code requires the Curve Fitting Toolbox in Matlab

%%  Equation 4: y = mean(y)
%  Data for 'eq_one_fit' fit:
%      X Input : x (trial number, 1:60)
%      Y Output: ERP amplitude of choice
%  Output:
%      fitresult : mean(y).
%      bic: Bayesian Information Criterion, used for model comparison


[xData, yData] = prepareCurveData( x, y );

coeffs = mean(y);
sse=(size(x,1)-1)*var(y);
bic = BIC_compute(sse, size(x,1), size(coeffs,2));

for t = 1:size(x,1)
fitresult(t) = mean(y);
end


function [fitresult, gof, bic] = eq_one_fit(x, y)

%% ERP modelling code relative to:
% Mancini F, Pepe A, Bernacchia, A, Di Stefano G, Mouraux A, Iannetti GD. (2017)
% Characterising the short-term habituation of event-related evoked potentials
% bioRxiv doi:10.1101/153387

% written in Matlab R2016b by F Mancini, fm456@cam.ac.uk
% This code requires the Curve Fitting Toolbox in Matlab


%%  Equation 1: y = a + b/x
%  Data for 'eq_one_fit' fit:
%      X Input : x (trial number, 1:60)
%      Y Output: ERP amplitude of choice
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%      bic: Bayesian Information Criterion, used for model comparison


[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a+b/x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.1 0.1];
opts.Lower = [0 0 ];
opts.Upper = [5 5 ];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
coeffs = coeffvalues(fitresult);
bic = BIC_compute(gof.sse, size(x,1), size(coeffs,2));



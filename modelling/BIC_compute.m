function bic = BIC_compute(SS,N,K)

% BIC, Schwarz Bayesian information criterion for curve-fitting model comparison
%   BIC(SS,N,K) where SS is the sum of squared residuals (as returned by e.g. lsqcurvefit)
%   N is the number of data points and K is the number of coefficients. 
%   Computes and returns the corrected AIC score for the model.
%   
%   The model with the lowest BIC value is the best fit.

  

K  = K + 1; % additional degree-of-freedom is model!

bic = (N .* log(SS./N)) + (K .* log(N));
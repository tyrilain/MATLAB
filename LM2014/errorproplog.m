%ERRORPROPLOG  Calculate the means and errors of F=log2(X)
%
%   [MEANS ERRORS] = ERRORPROPLOG(X,X_ERROR)
%       Where X is a vector of means and X_ERROR is a same length vector of
%       error estimates.
%
% Tara Martin, 2014-02-17


function [means errors] = errorproplog(X, X_error)

means = log2(X); % take base 2 log of means

% using the formula error_f = X_error/(ln(2)*X)
errors = X_error./(log(2)*X); % log(2) is the natural log of 2



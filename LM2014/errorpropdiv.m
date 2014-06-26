%ERRORPROPDIV  Calculate the mean and uncertainty of F=X/Y, where X and Y
%are same length vectors and we are really calculating f=x/y element by
%element.
%
%   [MEANS ERRORS] = ERRORPROPDIV(X,X_ERROR,Y,Y_ERROR)
%       Where X and Y are vectors of means and X_SE and Y_SE are the
%       corresponding error measurements.
%
% Tara Martin, 2014-02-17


function [means errors] = errorpropdiv(X, X_error, Y, Y_error)

% convert error measurements to relative error to make calculation easier
X_rel_error = X_error./X;
Y_rel_error = Y_error./Y;

% calculating the ratio of means is straight forward
means = X./Y;

% using the formula error_f = f*sqrt(rel_error_x^2 + rel_error_y^2)
errors = means.*sqrt(X_rel_error.^2+Y_rel_error.^2);



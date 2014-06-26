%ERRORPROP  Calculate the mean and uncertainty of F=X/Y, where X and Y
%are same length vectors and we are really calculating f=x/y element by
%element.
%
%   [MEAN ERROR] = ERRORPROP(FUNCTION,X,X_ERROR,Y,Y_ERROR)
%       Where FUNCTION is the mathematical function to perform, X and Y are 
%       vectors of means and X_SE and Y_SE are the corresponding error 
%       measurements.
%
%   FUNCTION options include 'add', 'subtract', 'multiply' and 'divide'.
%
% Tara Martin, 2014-02-17
% Edited: 2014-03-09    Added additional functions


function [means errors] = errorprop(func, X, X_error, Y, Y_error)

% calculate relative error to make multiplication or division easier
% using the formula error_f = f*sqrt(rel_error_x^2 + rel_error_y^2)
X_rel_error = X_error./X;
Y_rel_error = Y_error./Y;


switch lower(func)
    case 'add'
        means = X + Y;
        errors = sqrt(X_error.^2+Y_error.^2);
    case 'subtract'
        means = X - Y;
        errors = sqrt(X_error.^2+Y_error.^2);
    case 'multiply'
        means = X.*Y;
        errors = means.*sqrt(X_rel_error.^2+Y_rel_error.^2);
    case 'divide'
        means = X./Y;
        errors = means.*sqrt(X_rel_error.^2+Y_rel_error.^2);
    otherwise
        disp('Unknown function.')
        means = nan;
        errors = nan;
end



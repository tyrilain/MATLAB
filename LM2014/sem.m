function y = sem(x, w, dim)
%SEM Standard error of the mean.
%   For vectors, Y = SEM(X) returns the standard error.  For matrices,
%   Y is a row vector containing the standard error of each column.  For
%   N-D arrays, SEM operates along the first non-singleton dimension of X.
%
%   SEM normalizes Y by (N-1), where N is the sample size.  This is the
%   sqrt of an unbiased estimator of the variance of the population from
%   which X is drawn, as long as X consists of independent, identically
%   distributed samples.
%
%   Y = SEM(X,1) normalizes by N and produces the square root of the second
%   moment of the sample about its mean.  SEM(X,0) is the same as SEM(X).
%
%   Y = SEM(X,FLAG,DIM) takes the standard error along the dimension
%   DIM of X.  Pass in FLAG==0 to use the default normalization by N-1, or
%   1 to use N.
%
%   Example: If X = [4 -2 1
%                    9  5 7]
%     then sem(X,0,1) is [3.5355 4.9497 4.2426] and std(X,0,2) is [3.0
%                                                                  2.0]
%   Class support for input X:
%      float: double, single
%
%   See also STD, MEAN, VAR, MEDIAN, CORRCOEF.

%   Based on STD and VAR functions: Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 5.25.4.2 $  $Date: 2004/03/09 16:16:30 $

% Call var(x,flag,dim) with as many of those args as are present.

if nargin < 2 || isempty(w), w = 0; end
if nargin < 3
    % The output size for [] is a special case when DIM is not given.
    if isequal(x,[]), y = NaN(class(x)); return; end

    % Figure out which dimension sum will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

n = size(x,dim);
y = sqrt(var(x,w,dim)) ./ sqrt(n);
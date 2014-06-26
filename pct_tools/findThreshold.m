%FINDTHRESHOLD find the value of a distribution that is the mode + 1
%standard deviation
%
%    T = FINDTHRESHOLD(X)
%
%    X is a vector of data values
%    
% Zeba Wunderlich, 2010-07-26


function t = findThreshold(x)

x = x(:);
edges = linspace(min(x), max(x), length(x)/50);
[n, bin] = histc(x, edges);
if (mode(bin) < numel(edges))
    t = (edges(mode(bin)) + edges(mode(bin)+1))/2+std(x);
else
    t = (edges(mode(bin))+std(x));
end

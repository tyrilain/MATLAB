%HORIZONTAL_ERRORBAR   Compute plot data for horizontal error bars
% [xb,yb] = horizontal_errorbar(x,y,u,teesize)
% plot(xb,yb)

function [xb,yb] = horizontal_errorbar(x,y,u,tee)
% Copied from matlab/specgraph/errorbar.m

if nargin<4
   tee = (max(y(:))-min(y(:)))/100;  % make tee .02 y-distance for error bars
end

x = reshape(x,1,numel(x));
y = reshape(y,1,numel(y));
u = reshape(u,1,numel(u));
I = isnan(x) | isnan(y) | isnan(u);
if all(I)
   % nothing to draw
   xb = [];
   yb = [];
   return
end
x(I) = [];
y(I) = [];
u(I) = [];

n = numel(y);

yl = y - tee;
yr = y + tee;
xtop = x + u;
xbot = x - u;

yb = zeros(9,n);
yb(1,:) = y;
yb(2,:) = y;
yb(3,:) = NaN;
yb(4,:) = yl;
yb(5,:) = yr;
yb(6,:) = NaN;
yb(7,:) = yl;
yb(8,:) = yr;
yb(9,:) = NaN;

xb = zeros(9,n);
xb(1,:) = xtop;
xb(2,:) = xbot;
xb(3,:) = NaN;
xb(4,:) = xtop;
xb(5,:) = xtop;
xb(6,:) = NaN;
xb(7,:) = xbot;
xb(8,:) = xbot;
xb(9,:) = NaN;

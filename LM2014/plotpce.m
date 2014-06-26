%PLOTPCE   Display a pointcloud object
%   PLOTPCE(PC) displays the pointcloud object PC as a 2D plot in a/p-phi coordinates.
%   Expression values are mapped onto the scatter plot as different colors.
%
%   PLOTPCE(PC, AP_RANGE) displays the pointcloud object PC as a 2D plot in 
%   a/p-phi coordinates where AP_RANGE is of the form [ap_min ap_max].
%
%   PLOTPCE(PC, AP_RANGE, S) sets a custom cell marker size in pt^2. Useful
%   for adjusting marker size so cells don't overlap in figures.
%
% Tara Martin, 2014-03-25 (extracted from pointcloud/disp.m by Cris Luengo)
% Edits:
% 2014-04-10: added custom AP axis option.
%

function plotpce(pc, xrange, s)

x = pc.x; %pull out ap-axis values
coords = pc.apcoords;   %pull out ap coordinate system values
phi = atan2(coords(:,3),coords(:,2));   %convert from polar to linear coordinates

if nargin <2 || isempty(xrange)
    xrange = [min(x),max(x)];
end
if nargin <3
    s=10;
end

%plot pointcloud
scatter(x,phi,s,double(pc),'filled');

%set axis properties
set(gca,'xlim',xrange,'ylim',[-pi,pi],'ydir','reverse','ytick',[-pi 0 pi],...
         'yticklabel',{'D','V','D'},'dataaspectratio',[diff(xrange),1.2*2*pi,1]); 
     %edited aspect ratio to be not square, closer to actual embryo
set(gca,'box','on','visible','on');
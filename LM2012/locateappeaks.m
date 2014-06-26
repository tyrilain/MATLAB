%LOCATEAPEAKS   Determine the location of peaks of a/p stripes.
%   POS = LOCATEAPEAKS(VALUES,INIT) returns the a/p location
%   of stripe peaks in VALUES close to INIT. VALUES is an array (or
%   image) in which each row represents a resampled a/p strip, such as the
%   VALUES array returned by EXTRACTPATTERN. INIT is an array with
%   approximate locations, the returned locations in POS will be close
%   to those in INIT.
%
%   Locations are in egglength ratio (between 0 and 1, for anterior and
%   posterior respectively).
%
%   The rows in VALUES should be properly smoothed (eg using gaussf).
%
%   Example:
%   VALUES = 16-by-100 array
%   INIT = 16-by-4 array
%   POS = LOCATEAPEAKS(VALUES,INIT)

% Tara Martin, May 2012
% Based on LOCATEAPBOUNDARIES by
% Cris Luengo, somewhere in 2005

function [pos vals] = locateappeaks(onedim,init);

onedim = double(onedim);
N = size(onedim,1);        % number of a/p strips the embryo is divided in
L = size(onedim,2);        % number of samples along each a/p strip
M = size(init,2);          % number of peaks to measure stuff for
if size(init,1)~=N
   error('Wrong sizes for inputs');
end

pos = nan(N,M);
vals = nan(N,M);

extrema = double(onedim == dilation(onedim,3));
% extrema = double(onedim == dilation(onedim,[3,0]));
% extrema = extrema - double(onedim == erosion(onedim,[3,0]));
xx = (0:(L-1))/(L-1);

for ii=1:N       % ii is the current d/v strip
    f = onedim(ii,:);  % single trace
    maxi = find(extrema(ii,:)==1);   % Index of local maxes
    % filter for peaks above background (ie mean)
%     maxi = maxi(f(maxi)>mean(f));
    maxi = maxi(f(maxi)>0.1);
    e = xx(maxi); % get ap positions of these
    
    for jj=1:M    % jj is the current expression stripe
      [tmp,I] = min(abs(e-init(ii,jj)));
      E = maxi(I);  % this is the index of the nearest peak
      % finding peak by center of mass of top 10% of stripe
      v = f(E);
      l = double(label(f>0.9*v,1));
      I = find(l==l(E));
      pos(ii,jj) = sum(f(I).*xx(I))/sum(f(I));
      vals(ii,jj) = v;
   end
end

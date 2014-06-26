function cmap = makecolmap(mapnum,n);
if nargin<2
   n = size(get(gcf,'colormap'),1);
end
      co0 = [ 0.0, 0.0, 0.0 ];
switch mapnum
   case 1 % blue
      co1 = [ 0.0, 0.0, 1.0 ];
      co2 = [ 0.0, 0.8, 1.0 ];
   case 2 % red
      co1 = [ 1.0, 0.0, 0.0 ];
      co2 = [ 1.0, 0.8, 0.0 ];
   case 3 % green
      co1 = [ 0.0, 0.8, 0.0 ];
      co2 = [ 0.7, 0.8, 0.0 ];
   otherwise
      error('I can''t make more than 3 different color maps!')
end
if n > 3
   %m = ceil(n*0.7);
   %if mapnum==1
      m = ceil(n*0.6); % More emphasis on second half of color map.
   %end
   cmap1 = [linspace(co0(1),co1(1),m)',    linspace(co0(2),co1(2),m)',    linspace(co0(3),co1(3),m)'];
   cmap2 = [linspace(co1(1),co2(1),n-m+1)',linspace(co1(2),co2(2),n-m+1)',linspace(co1(3),co2(3),n-m+1)'];
   cmap = [cmap1;cmap2(2:end,:)];
else
   cmap = [co0;co1;co2];
   cmap = cmap(1:n,:);
end

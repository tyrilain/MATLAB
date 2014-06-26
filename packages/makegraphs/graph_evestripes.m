% Plot eve stripes averaged per cohort

clear all
close all

% Settings that make this script specific to EVE:
dbname = 'eve_stripes_data';
imgname = 'eve_stripes';

% If WIDEGRAPH is TRUE, we're folding left and right sides.
widegraph = false;
% widegraph = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What's below should be the same for any A/P patterning gene:

load(dbname)
% results == {filename,cohort,init,pos}

cohorts  = unique([results{:,2}]);
Ncohorts = length(cohorts);
Npoints  = size(results{1,4},1);
Nlines   = size(results{1,4},2);

% Compute averages for stripe locations and intensities
if widegraph
   NPsz = Npoints;
else
   NPsz = Npoints/2+1;
end
location = [];
location.pos = nan(NPsz,Nlines,Ncohorts); % average coordinates in %EL
location.std = nan(NPsz,Nlines,Ncohorts); %    standard deviations
location.err = nan(NPsz,Nlines,Ncohorts); %    95% confidence intervals
Nembryos = zeros(Ncohorts,1);
legendstr = cell(Ncohorts,1);
for ii=1:Ncohorts
   I = find([results{:,2}] == cohorts(ii));
   N = length(I);
   Nembryos(ii) = N;
   if N>0
      loc = zeros(Npoints,Nlines,N);
      for kk=1:N
         loc(:,:,kk) = results{I(kk),4};
      end
      if ~widegraph
         loc = folddata(loc);
      end
      [m,s,e] = averages(loc);
      location.pos(:,:,ii) = m;
      location.std(:,:,ii) = s;
      location.err(:,:,ii) = e;
   end
   legendstr{ii} = ['cohort #',num2str(cohorts(ii))];
end
Npoints  = NPsz;

% Some more settings

linewidth = 1;
fontsize1 = 14;
fontsize2 = 20;

% Do the plotting
l  = permute(location.pos,[2,1,3]);
l2 = permute(location.err,[2,1,3]);
if widegraph
   stripe_graph(l,l2,legendstr);
else
   stripe_graph(l,l2,legendstr,'half');
end
set(gca,'linewidth',linewidth,'fontsize',fontsize1);
set(findobj(gca,'type','line'),'linewidth',linewidth);
set(get(gca,'XLabel'),'fontsize',fontsize1)
set(get(gca,'YLabel'),'fontsize',fontsize1)
print('-depsc2','-r300','-loose',[imgname,'.eps'])
system(['convert -size 200x200 ',imgname,'.eps',' -resize 200x200 ',imgname,'_thumb.png']);

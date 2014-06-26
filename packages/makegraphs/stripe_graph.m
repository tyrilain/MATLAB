function hout = stripe_graph(m_lines,u_lines,legendstr,mode,highlight,alpha)
% Plotting AP stripes in an unrolled view, with errorbars and a legend
%
% STRIPE_GRAPH(M_LINES,U_LINES,LEGENDSTR,MODE,HIGHLIGHT,ALPHA)
%
% M_LINES is an Nlines x Npoints x Ncohorts array. M_LINES(LINE,:,COHORT)
% gives the mean location of Npoints around stripe number LINE for cohort
% COHORT. U_LINES gives the corresponding extent of the error bars.
%
% LEGENDSTR is a cell array with Ncohorts strings.
%
% MODE is a string or string cell array with options:
%   - 'sparse' means only 1/3 of the data points will have errorbars. By
%     default, each point will have errorbars.
%   - 'micron' means the x-axis is in micron. It will range from -200 to
%     200. By default it is in %EL, from 0 to 1.
%   - 'fixed' means the y-axis will not be re-arranged. Use this if your
%     data already has the ventral in the middle and dorsal point duplicated.
%   - 'half' means the data goes from ventral to dorsal.
%   - combine 'half' and 'fixed' if your data goes from dorsal to ventral.
%   - 'alternatecolors' makes the background black for a different contrast.
%
% HIGHLIGHT is a cell array with 2 x Npoints x Ncohorts arrays. Each
% of these indicate a band that will be colored. HIGHLIGHT{i} can also have
% 2 x Npoints x Ncohorts+1 elements. In that case, HIGHLIGHT{i}(:,:,END)
% indicates a band to be drawn in grey (it will be behind the other bands).
%
% ALPHA is the alpha value for the highlight. Set to 1 to not use
% transparancy. Set to 0 (the default) to use hatching instead of solid
% coloring.


% check which inputs are specified, and select mode
if nargin<6
   alpha = 0;
end
if nargin<5
   highlight = {};
end
if nargin<4
   mode = {};
else
   if ischar(mode)
      mode = {mode};
   elseif ~iscellstr(mode)
      error('MODE should be a string or a cell array with strings.')
   end
end
if any(strcmp(mode,'half'))
   halfembryo = true;
else
   halfembryo = false;
end
if any(strcmp(mode,'alternatecolors'))
   alternatecolors = true;
else
   alternatecolors = false;
end
if nargin<3
   legendstr = {};
end
if nargin<2
   u_lines = zeros(size(m_lines));
end
if ~isequal(size(m_lines),size(u_lines))
   error('M_LINES and U_LINES not of same size')
end
Nlines = size(m_lines,1);
Npoints = size(m_lines,2);
Ncohorts = size(m_lines,3);

for ii=1:length(highlight)
   [sz1,sz2,sz3] = size(highlight{ii});
   if sz1~=2 || sz2~=Npoints || ( sz3~=Ncohorts && sz3~=Ncohorts+1 )
      error('HIGHLIGHT array is of incorrect size')
   end
end

if ~any(strcmp(mode,'fixed'))
   if halfembryo
      % Turn the data so that the dorsal side is up
      m_lines = m_lines(:,[end:-1:1],:);
      u_lines = u_lines(:,[end:-1:1],:);
      for ii=1:length(highlight)
         highlight{ii} = highlight{ii}(:,[end:-1:1],:);
      end
   else
      % Turn the data so that the embryo is cut on the dorsal side rather
      % than the ventral side
      cutpoint = floor(Npoints/2)+1;
      m_lines = m_lines(:,[cutpoint:Npoints,1:cutpoint],:);
      u_lines = u_lines(:,[cutpoint:Npoints,1:cutpoint],:);
      for ii=1:length(highlight)
         highlight{ii} = highlight{ii}(:,[cutpoint:Npoints,1:cutpoint],:);
      end
      Npoints = Npoints+1;
   end
end

% set line colors for each cohort
if Ncohorts==2
   cols = [0,0,1;1,0,0];
elseif Ncohorts==3
   % Like N=2, third color is green.
   cols = [0,0,1;1,0,0;0,0.8,0];
elseif Ncohorts==4
   % 4-color JET
   cols = [0,0,1;0,0.9,0.5;0.8,0.7,0;1,0,0];
else
   cols = jet(200);
   cols = cols(round(linspace(25,175,Ncohorts)),:);
end

% selecting hatching styles
hatchangle = [];
while Ncohorts>length(hatchangle)
   hatchangle = [hatchangle,45,135,0,90,45/2,5*45/2,3*45/2,7*45/2];
end

% make a figure with appropriate graphic settings and legends
figure, hold on, box on, set(gcf,'InvertHardcopy','off','color','white');
if halfembryo
   set(gcf,'position',[500,280,800,500],'paperpositionmode','auto')
else
   set(gcf,'position',[500,280,800,800],'paperpositionmode','auto')
end
if alternatecolors
   set(gca,'color','black')
end
if any(strcmp(mode,'micron'))
   set(gca,'xlim',[-200,200])
   xlabel('AP position (micron from center of mass)');
else
   set(gca,'xlim',[0,1],'xtick',0:0.1:1,'xticklabel',0:10:100)
   xlabel('AP position (% egglength)')
end
if halfembryo
   segments = 1:(Npoints-1)/2:Npoints;
   segmentnames = {'dorsal','side','ventral'};
else
   segments = 1:(Npoints-1)/4:Npoints;
   segmentnames = {'dorsal','left','ventral','right','dorsal'};
end
set(gca,'ylim',[0.5,segments(end)+0.5],'ytick',segments,'yticklabel',segmentnames,'ydir','reverse')
ylabel('DV location')

if alternatecolors
   c = [0.4,0.3,0.4];
else
   c = [1,1,1]*0.9;
end

% draw the highlights
for jj=1:length(highlight)
   hl = highlight{jj};
   if size(hl,3) > Ncohorts
      % This band should be in grey behind all the other stuff
      ii = Ncohorts+1;
      I = isnan(hl(1,:,ii)) | isnan(hl(2,:,ii));
      J = [find([~I(1),I(1:end-1) & ~I(2:end)]),length(I)+1];
      for kk=1:length(J)-1
         x = [hl(1,J(kk):J(kk+1)-1,ii),fliplr(hl(2,J(kk):J(kk+1)-1,ii))];
         y = [J(kk):J(kk+1)-1,fliplr(J(kk):J(kk+1)-1)];
         K = [I(J(kk):J(kk+1)-1),fliplr(I(J(kk):J(kk+1)-1))];
         x(K) = [];
         y(K) = [];
         patch(x,y,c,'EdgeColor','none');
      end
   end
end
for jj=1:length(highlight)
   hl = highlight{jj};
   for ii=1:Ncohorts
      I = isnan(hl(1,:,ii)) | isnan(hl(2,:,ii));
      J = [find([~I(1),I(1:end-1) & ~I(2:end)]),length(I)+1];
      for kk=1:length(J)-1
         x = [hl(1,J(kk):J(kk+1)-1,ii),fliplr(hl(2,J(kk):J(kk+1)-1,ii))];
         y = [J(kk):J(kk+1)-1,fliplr(J(kk):J(kk+1)-1)];
         K = [I(J(kk):J(kk+1)-1),fliplr(I(J(kk):J(kk+1)-1))];
         x(K) = [];
         y(K) = [];
         c = cols(ii,:);
         if alpha>0
            if alternatecolors
               c = 1-(1-c)/2;
            else
               c = 1-(1-c)/6;
            end
            patch(x,y,c,'FaceAlpha',alpha,'EdgeColor','none');
         else
            if ~alternatecolors
               c = 1-(1-c)/2;
            end
            h = patch(x,y,'k','EdgeColor','none','FaceColor','none');
            hatch(h,hatchangle(mod(ii-1,length(hatchangle))+1),c,'-',7,0.5)
         end
      end
   end
end

% this, um, doesn't do anything?
if alternatecolors
   %cols = 1-(1-cols)/6;
end

% finally, draw our stripe boundaries and error bars!
h = zeros(Ncohorts,1);
for ii=1:Ncohorts
   x = m_lines(1,:,ii);
   if any(isfinite(x))
      y = 1:Npoints;
      h(ii) = plot(x,y,'color',cols(ii,:),'marker','.');
      if any(strcmp(mode,'sparse'))
         I = 1:3:Npoints;
         x = x(I);
         y = y(I);
         u = u_lines(1,I,ii);
      else
         u = u_lines(1,:,ii);
      end
      [x,y] = horizontal_errorbar(x,y,u,0.15);
      plot(x,y,'color',cols(ii,:));
   else
      h(ii) = plot(nan,nan,'color',cols(ii,:),'marker','.');
   end
   for jj=2:Nlines
      x = m_lines(jj,:,ii);
      if any(isfinite(x))
         y = 1:Npoints;
         plot(x,y,'color',cols(ii,:),'marker','.');
         if any(strcmp(mode,'sparse'))
            I = mod(jj-1,3)+1:3:Npoints;
            x = x(I);
            y = y(I);
            u = u_lines(jj,I,ii);
         else
            u = u_lines(jj,:,ii);
         end
         [x,y] = horizontal_errorbar(x,y,u,0.15);
         plot(x,y,'color',cols(ii,:));
      end
   end
end

% add a legend
if halfembryo
   legend(h,legendstr,'Location','NorthOutside','color','white')
else
   legend(h,legendstr,'Location','SouthWest','color','white')
end

% optionally return figure handle
if nargout>0
   hout = h;
end

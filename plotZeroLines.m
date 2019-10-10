function plotZeroLines( subp_y, subp_x, varargin )
% plotZeroLines( subplots_y, subplots_x, ... )
% plot lines at x==0, y==0
% sgm 2018

if exist('subp_y','var') % subplots?  do it for each one.
  for i=1:subp_x
    for j=1:subp_y
      subplot( subp_y,subp_x, i + (subp_x)*(j-1) )
      plotlines();
    end
  end
else  % no subplots, just do it once
  plotlines()
end    


function plotlines
  held = ishold;
  hold on;
  plot(xlim,[0 0],':','linewidth',2);
  plot([0 0],ylim,':','linewidth',2);
  xl = xlim;
  yl = ylim;
  if xl==yl
    plot(xl,xl,':','linewidth',2)
  end
  if ~held
    hold off;
  end


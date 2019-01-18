function varargout = imagep(p, ylabels, xlabels,varargin)
%   handle =  imagep(p)
%   handle =  imagep(p, ylabels, xlabels)
%  plot a matrix of p-values, e.g. from a correlation.
% 
%  calls imagesc(p)
%  but rescales colours for p-value interpretation.
%  so the colour is -log(p), and the colorbar shows a log-scale.
%  
%  will accept signed p-value (positive or negative) if a signed contrast was performed.
%  default alpha = 0.05 - edit code to change this!
%
% sgm 2012
ALPHA  = 0.05;
ROTATE = true; % make x axis labels run vertically?
sign=(p>0)-(p<0); p=abs(p); % keep track of sign
sign(sign==0)=1; % p-values of zero are taken to be +0.00001
if any(p(:)< ALPHA), p(p>ALPHA)=1; end
if 0
  p(p<eps)=eps; % log(0) = inf, so use log(eps) instead -- gives approximately 16.
else
  p(p<eps)=1; % count zeros as p=1
end
imval = -log(p)/log(10);
imval = imval .* sign;
if nargout>0
  varargout{1:nargout} = imagesc(imval);
else
  imagesc(imval);
end
if exist('xlabels','var') && ~isempty(xlabels)
  if ROTATE
    set(gca,'XAxisLocation','top');
  end
  set(gca,'xtick',1:length(xlabels),'xticklabel',xlabels);
  if ROTATE
    rotateXLabels(gca,90);
    set(gca,'XAxisLocation','bottom');
  end
end
if exist('ylabels','var') && ~isempty(ylabels)
  set(gca,'ytick',1:length(ylabels),'yticklabel',ylabels);
end
if any(sign(:)<0)
  caxis([-2 2]); 
else
  caxis([0 2]); 
end
h=colorbar; 
set(h,'yticklabel',arrayfun(@(x) ...
  sprintf('%0.2g',x), 10.^-abs(get(h,'ytick')) ...
  , 'uniformoutput',0));

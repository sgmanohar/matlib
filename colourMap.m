function [ col ] = colourMap( index, maximum, varargin )
% colourMap(index, maximum) 
%   an easy and safe way to index into the colour palette.
%   Returns a colour from the current colour map. 
%   'index' should lie from 1 to  maximum.
% colourMap with no arguments:
%   sets the current figure's default line colours using the current
%   colormap, i.e. 
%   set(gcf,'DefaultAxesColorOrder', colormap)
% colourMap(maximum)
%   set the current figure's default line colours using a random palette
%   selected from 'othercolor', from 1 to maximum. 
%   Also set the colormap.
% colourMap( palette )
%   set the line colours from the given palette
if ~exist('index','var')  % no arguments? set axes color order to color map
  set(gcf, 'DefaultAxesColorOrder', colormap);
  set(gca, 'ColorOrder', colormap);
  return;
end

if exist('maximum','var') && ~isnumeric(maximum) 
  varargin=[varargin {maximum}]; 
  clear maximum; 
end

i=find(strcmpi('emphasize',varargin));
if length(i)==1
  emphasize = varargin{i};
else emphasize = [];
end

if ~exist('maximum','var') % one argument:
  if isscalar(index) % scalar? create a default colormap with n lines
    if ~exist('othercolor','file') || 0 % Always use Jet palette?
      col = jet(index);
    else % use a random 'othercolor' palette
      col = othercolor( 1+floor(rand*404), index);
    end
    set(gcf, 'DefaultAxesColorOrder', col );
    set(gca, 'ColorOrder', col );
    colormap(col);
  elseif ismatrix(index) && size(index,2)==3 % is it a colourmap? just use that
    set(gcf, 'DefaultAxesColorOrder', index );
    set(gca, 'ColorOrder', index );
  elseif ischar(index)
    switch index
      case 'centred'
        col = centred_red_blue();
      otherwise
        error('unrecognised palette %s', index);
    end
    set(gcf, 'DefaultAxesColorOrder', col );
    set(gca, 'ColorOrder', col );
    colormap(col);
  end
else % two arguments: index and maximum for current colormap
  c=colormap;
  index=max(min(index,maximum),1); % force index into range 1:maximum
  ix = floor((size(colormap,1)-1) * (index-1)/maximum) +1; % index into colour palette
  col=c(ix,:);
end

if ~isempty(emphasize)
  
end

return




function Y = centred_red_blue
% make it symmetrical about zero
a=caxis;
mx = max(-a(1), a(2));
caxis([-mx, mx]);

Y = [
  
         0         0    1.0000
         0    0.0625    1.0000
         0    0.1250    1.0000
         0    0.1875    1.0000
         0    0.2500    1.0000
         0    0.3125    1.0000
         0    0.3750    1.0000
         0    0.4375    1.0000
         0    0.5000    1.0000
         0    0.5625    1.0000
         0    0.6250    1.0000
         0    0.6875    1.0000
         0    0.7500    1.0000
         0    0.8125    1.0000
         0    0.8750    1.0000
         0    0.9375    1.0000
         0    1.0000    1.0000
    0.0695    0.9980    0.9980
    0.1389    0.9961    0.9961
    0.2084    0.9941    0.9941
    0.2779    0.9922    0.9922
    0.3473    0.9902    0.9902
    0.4168    0.9882    0.9882
    0.4863    0.9863    0.9863
    0.5557    0.9843    0.9843
    0.6252    0.9824    0.9824
    0.6947    0.9804    0.9804
    0.7641    0.9784    0.9784
    0.8336    0.9765    0.9765
    0.9031    0.9745    0.9745
    0.9725    0.9725    0.9725
    0.9743    0.9743    0.9118
    0.9760    0.9760    0.8510
    0.9777    0.9777    0.7902
    0.9794    0.9794    0.7294
    0.9811    0.9811    0.6686
    0.9828    0.9828    0.6078
    0.9846    0.9846    0.5471
    0.9863    0.9863    0.4863
    0.9880    0.9880    0.4255
    0.9897    0.9897    0.3647
    0.9914    0.9914    0.3039
    0.9931    0.9931    0.2431
    0.9949    0.9949    0.1824
    0.9966    0.9966    0.1216
    0.9983    0.9983    0.0608
    1.0000    1.0000         0
    1.0000    0.9412         0
    1.0000    0.8824         0
    1.0000    0.8235         0
    1.0000    0.7647         0
    1.0000    0.7059         0
    1.0000    0.6471         0
    1.0000    0.5882         0
    1.0000    0.5294         0
    1.0000    0.4706         0
    1.0000    0.4118         0
    1.0000    0.3529         0
    1.0000    0.2941         0
    1.0000    0.2353         0
    1.0000    0.1765         0
    1.0000    0.1176         0
    1.0000    0.0588         0
    1.0000         0         0
];

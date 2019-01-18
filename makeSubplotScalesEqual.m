function makeSubplotScalesEqual(h,w,IDX,  FIG )
%  makeSubplotScalesEqual(h,w,IDX, [FIG])
% 
%   go through each of the subplots given in IDX,
%   calculate the smallest scale that accommodates them all,
%   and rescale them all to the this size.
%   if IDX is omitted, use all the subplots (i.e. 1:h*w)
% 
%   if FIG is specified, it is a vector the same size as IDX,
%   and the specified subplots from a different figure are used.
%   note: the other figure must have the same shape subplots!
%
% SGM 2011
if ~exist('IDX','var')
    IDX=1:(h*w);
end
if isempty(IDX), return; end

% read all axes limits
for(i=1:length(IDX))
  if(exist('FIG','var')) figure(FIG(i)); end
  subplot(h,w,IDX(i));
  xl(i,:)=xlim();
  yl(i,:)=ylim();
  zl(i,:)=zlim();
  cl(i,:)=caxis();
end
% find largest range that encapsulates all of them
nxl = [ min(xl(:,1)) max(xl(:,2)) ]; % new x limits
nyl = [ min(yl(:,1)) max(yl(:,2)) ]; % new y limits
nzl = [ min(zl(:,1)) max(zl(:,2)) ]; % new z limits
ncl = [ min(cl(:,1)) max(cl(:,2)) ]; % new colour limits

% set new axes limits
for(i=1:length(IDX))
  if(exist('FIG','var')) figure(FIG(i)); end
  subplot(h,w,IDX(i))
  xlim(nxl);
  ylim(nyl);
  zlim(nzl);
  caxis(ncl);
end



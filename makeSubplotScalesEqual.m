function makeSubplotScalesEqual(h,w,IDX,  FIG, varargin )
%  makeSubplotScalesEqual(h,w,IDX, [FIG])
% 
%   go through each of the subplots given in IDX,
%   calculate the smallest scale that accommodates them all,
%   and rescale them all to the this size.
%   This will apply to: xlim, ylim, zlim, and caxis.
%   
%   if IDX is omitted, use all the subplots (i.e. 1:h*w)
%   e.g.  makeSubplotScalesEqual( 2,2 );  % make the scales of all 4 subplots equal 
% 
%   if FIG is specified, it is a vector the same size as IDX,
%   and the specified subplots from a different figure are used.
%   note: the other figure must have the same shape subplots!
%   e.g. makeSubplotScalesEqual( 2,2, [1 2 1 2] , [1 1 2 2] ); 
%    makes the scales of the upper two subplots on both figures equal.
%
%   if FIG is 'across', then the next parameter must be a list of figure
%   indices. Matching subplots across those figures are made equal. 
%   e.g. makeSubplotScalesEqual( 2,2, 1:4, 'across', 1:2 )
%    make the scales of subplot 1 identical across the two figures,  then
%    make the scales of subplot 2 identical across the two figures. 
% 
% SGM 2011-2019

if ~exist('IDX','var') || isempty(IDX) % default to all subplots
    IDX=1:(h*w);
end
if isempty(IDX), return; end % no subplots? maybe h or w are zero. 

% are we doing each subplot independently, across multiple figures? 
if exist('FIG','var') && strcmp(FIG, 'across') 
  FIG = varargin{1};      % get list of figures
  varargin=varargin(2:end);
  for i=1:length(IDX)     % call the function for each subplot separately. 
    makeSubplotScalesEqual(h,w,repmat(IDX(i),size(FIG)), FIG);
  end
  return
end

i=find(strcmpi('subplot',varargin),1);
if any(i)
  SUBPLOT = varargin{i+1};
  varargin(i:i+1)=[]; 
else
  SUBPLOT = @subplot;
end

% read all axes limits 
for(i=1:length(IDX))
  if(exist('FIG','var') && ~isempty(FIG) ) figure(FIG(i)); end
  SUBPLOT(h,w,IDX(i));
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
  if(exist('FIG','var') && ~isempty(FIG) ) figure(FIG(i)); end
  SUBPLOT(h,w,IDX(i))
  xlim(nxl);
  ylim(nyl);
  zlim(nzl);
  caxis(ncl);
end



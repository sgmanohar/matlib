function [resid betas] = mainSequence(info, DIM, varargin)
% [resid,  betas] = mainSequence( info, [DIM] )
% calculate residual peak velocity from main sequence, given an 'info'
% structure. if DIM is included, compute separately for each level of DIM.
% 
if ~exist('DIM','var'), 
  DIM = length(size(info.sRT))+1; % use a blank dimension
end
i = find(strcmpi(varargin,'amprange'));
if i>0
  amprange = varargin{i+1};
  varargin([i i+1])=[];
else
  amprange = [-inf inf];
end
i = find(strcmpi(varargin,'velrange'));
if i>0
  velrange = varargin{i+1};
  varargin([i i+1])=[];
else
  velrange = [-inf inf];
end


NS = size( info.sRT, DIM );
index = repmat( {':'}, max(DIM, ndims(info.sRT)) , 1 ); % vector of indices
sz = size( info.sRT ); 
sz(DIM) = 1; % create a size structure for the individual subject data
for i=1:NS
  index{DIM} = i;
  amp = info.sAmpl(index{:});
  vel = info.sSpd(index{:});
  bad = isnan(amp) | isnan(vel) ... % select non-nan trials
      | amp < amprange(1)  | amp > amprange(2) ... % out of range?
      | vel < velrange(1)  | vel > velrange(2) ; 
  [b,~,res] = regress( vel(~bad), [ amp(~bad) ones(size(amp(~bad))) ] );
  resid{i} = nan(sz);
  resid{i}(~bad) = res;
  betas(i,:) = b;
end

resid = nancat([2,DIM],resid);
resid = resid{1};

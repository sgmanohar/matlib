function  [Y varargout]=alignRight(X, varargin)
% Y = alignRight(X)
%   takes a matrix that is padded with nans at the right edge
%   and align so that the nans are on the left.
% [Y A0 B0] = alignRight(X, A,B,...)
%   matrices A and B are of same size as X. 
%   shift these matrices exactly the same amount as X, discarding any
%   elements to the far right of those matrices

if nargin>1
  for j=1:nargin-1
    if any(size(X)~=size(varargin{j}))
      error('all args must be same size')
    end
  end
end

w=size(X,2); % width w
for i=1:size(X,1) % for each row
  r=X(i,:);
  ix=find(~isnan(r),1,'last') ; % find last non-nan value in the row
  if isempty(ix), ix=0; end
  Y(i,:)= [nan(1,w-ix)  r(1:ix)]; % nan pad left.
  for j=1:nargin-1
    O{j}(i,:)=[nan(1,w-ix)  varargin{j}(i,1:ix)];
  end
end
if exist('O','var')
  varargout=O;
end
function [h,X]=surfsmooth(varargin)
% function surfsmooth(varargin)
% calls surf, but applies 20% smoothing and makes the surface look nice.
% 'smooth' - the size of smoothing window to use, as a fraction of the grid
%            size. Default 0.2;
% sgm 
X=varargin{1};
KERNEL = 'GAUSS';
SMOOTH = 0.2; 
AXIS3D = true; % make the axis isometric on all 3 dimensions?

i=find(strcmpi(varargin,'smooth')); if ~isempty(i)
  SMOOTH=varargin{i+1}; varargin([i i+1])=[]; 
end

if SMOOTH && ismatrix(X) && numel(X)>30 
  sz=min(size(X));sz=floor(sz*SMOOTH); % smoothing: one fifth of size
  if sz<2, warning('surface is very small, smoothing wont work'); end
  switch KERNEL
    case 'GAUSS' % gaussian kernel
      kernel = exp( - bsxfun(@plus, linspace(-3,3,sz+1).^2, linspace(-3,3,sz+1)'.^2 )  ) ;
  end
  kernel = kernel / sum(kernel(:)); % normalise to 1
  % nanconv calls conv2 under the hood. 
  X=nanconv( X, kernel, 'edge' );
end
varargin{1}=X;
h=surf(varargin{:},'EdgeColor','None');
if AXIS3D
  axis vis3d
end
camlight
lighting phong
shading interp


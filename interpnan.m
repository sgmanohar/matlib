function X = interpnan(DIM, X, varargin)
%   interpnan ( DIM, X ... )
% interpolate across nanvalues, along dimenson DIM.
% it finds series of nans, and calls interp1, assuming equal
% spacing of the values. The result is filled in place of 
% the nans.
% varargin is passed to interp1.
% alternative: 
%   interpnan ( DIM, X, V, Xq, ... )
% resamples the variables as well as interpolating nans.
% sgm 2016



if  nargin>3 && isnumeric(varargin{1}) && isnumeric(varargin{2})
  C=X; X=varargin{1}; Cq = varargin{2};   % treat the inputs as ( DIM, X, V, Xq )
  varargin(1:2) = [];
  useQ = true;
else useQ = false; 
end
sz=size(X);
X = shiftdim(X,DIM-1);
X = X(:,:);
failed_count = 0;
for i=1:size(X,2)
  notnanix = find(~isnan(X(:,i))); 
  nanix    = find( isnan(X(:,i)));
  if length(notnanix) > 2 % if there are at lesat 3 points to interpolate between
    % has the user provided an X and Xq parameter (the x-coordinates, called C here)?
    if useQ  % use the provided coordinates to interpolate
      X(:, i) = interp1( C(notnanix), X(notnanix, i), Cq, varargin{:} );
    else % use the linear index instead
      X(nanix, i) = interp1( notnanix, X(notnanix, i), nanix, varargin{:});
    end
  else % otherwise, nan the whole vector
    x(notnanix,i) = nan;
    failed_count = failed_count + 1;
  end
end
if failed_count
  warning('%g vectors had too few datapoints and were nanned',failed_count)
end
X = shiftdim(X, -(DIM-1));
X = reshape(X,sz); 

function [P, t] = poincare_map(X, plane)
% P = poincare_map( X [, plane] )
% given N-dimensional time-series data X,  
% find points of the time series that cross a given plane. 
%   X ( t, variables ) is a T x N matrix of N-dimensional states evolving over time
%   plane.norm = N-dim normal vector  (default [1,0,0,0...]
%   plane.dist = distance from origin (default 0)
% sanjay g manohar 2019

FAST_MODE = false; % calculate in matrix form, but do not interpolate points?

% define the poincare plane
if ~exist('plane','var')
  plane.norm = zeros(size(X,2),1);  % column vec
  plane.norm(1) = 1;
  plane.dist = 0;
end
plane.norm = plane.norm/norm(plane.norm); % ensure length = 1

% distance of each state point from the plane
% project the state onto the plane's normal vector (dot product)
projections = X * plane.norm - plane.dist;
if FAST_MODE 
  % select points where the plane is crossed
  crossings = projections(1:end-1) < 0 & projections(2:end) > 0 ;
  P = X(crossings,:); 
  if nargout>1
    t=find(crossings); 
  end
else % SLOW MODE: estimate each plane-crossing using linear interpolation
  P = [] ; % keep a list of plane-crossing points
  t = [] ; % keep a list of the times of crossing
  for i=2:size(X,1) % for each timepoint
    p1 = projections(i-1); % previous state's distance from plane
    p2 = projections(i);   % current state's distance from plane
    % is this a plane-crossing?
    if p2>0 && p1<0
      % linear interpolation to mix the two points, to estimate the crossing
      % point
      mix = -p1/(p2-p1); % how much of p2 to mix into p1?
      Xp = X(i-1,:) + mix * (X(i,:)-X(i-1,:)); % estimated plane crossing
      % add a row to P
      P = [ P;  Xp ];
      t = [ t;  i  ];
    end
  end
end
% note that, if you have some directions within the plane, you can 
% collapse P into a set of (N-1)-dimensional points
return
%% Example
% test it on rossler equation:
x=[ 2 0 0 ]; c = 4.16; a = 0.2;   % initial conditions and parameters
for i=1:10000;                    % for each timepoint
  dx = [ -x(2)-x(3), x(1)+a*x(2), a+x(3)*(x(1)-c) ]; % calculate rossler
  x=x+dx*0.02;    xx(i,:)=x;      % update and store
end; % then plot trajectory
subplot(3,3,1); plot(xx(:,1)); title 'X'; xlabel t
subplot(3,3,2); plot(xx(:,2)); title 'Y'; xlabel t
subplot(3,3,3); plot(xx(:,3)); title 'Z'; xlabel t
subplot(3,3,4); plot3(xx(:,1),xx(:,2),xx(:,3)); xlabel X; ylabel Y; zlabel Z;
P = poincare_map( xx );           % compute X-plane crossings
subplot(3,3,5); plot(P(:,2),P(:,3),'o-'); xlabel Y,ylabel Z;
subplot(3,3,6); plot( P(1:end-1,2), P(2:end,2) , 'o'); xlabel Y_{t-1}; ylabel Y_t; title 'poincare map'
subplot(3,3,7); plot(P(:,1),'o-'); title X; xlabel cycle; 
subplot(3,3,8); plot(P(:,2),'o-'); title Y; xlabel cycle
subplot(3,3,9); plot(P(:,3),'o-'); title Z; xlabel cycle


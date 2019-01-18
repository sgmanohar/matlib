function p = logpolar_gaussian( r0, t0, rs, ts, N )
% p = logpolar_gaussian( r0, t0, rs, ts, N )
%
% probability density over a 2d plane, in a log-polar gaussian distribution
% centred on radius = r0, angle theta = t0.
% example:   imagesc( logpolar_gaussian( 0.5, pi/4 ) )
%            Plots a gaussan pdf in log-polar coordinates, centred at 
%            radius = 0.5, angle = 45 deg
%
% The log polar gaussian is a common model for the shape of visual
% receptive fields in retinotopic coordinates.
% 
% the plane extends from -1 to +1, and N indicates the number of divisions
% of the x and y axes. rs and ts indicate standard deviation for the r and 
% theta coordinates. Default N=200, ts=0.1, rs=0.2.
% 
% sgm 2016
if ~exist('N','var'),   N  = 200; end   % use default values if not specified
if ~exist('ts','var'),  ts = 0.1; end
if ~exist('rs','var'),  rs = 0.2; end
if ~exist('t0','var'),  t0 = 0; end
if ~exist('r0','var'),  r0 = 0; end

xx = linspace(-1,1,200);  % create a grid of equally spaced x and y values
[x,y] = meshgrid(xx,xx); 
z = x+1j*y;               % complex number
t = angle(z);             % to convert to polar coords
r = abs(z); 
p = normpdf( sqrt( log(r/r0).^2 / rs^2  +  ... % compute probability for each point
                      (t-t0).^2 / ts^2     ) );
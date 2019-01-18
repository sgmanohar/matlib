function z = gabor_matrix( gridSize, spatialFreq, orientation, falloff, phase ) 
% z = gabor_matrix( gridSize, spatialFreq, orientation ) 
% gridSize    - size of the output array, e.g. 100 gives a 100x100 matrix.
% spatialFreq - cycles per gridsize. 1 means one cycle across the whole
%               grid. 
% orientation - angle of the gabor in radians. 0 means horizontal; pi/2 is vertical. 
% falloff     - distance units (where size of grid = 1) where amplitude 
%               should fall off by 1 exponential unit.
%               A value of 1 means the amplitude will be e^-1 = 0.37 at the
%               edge of the grid.
% phase       - phase of the cosine wave, at the grid's centre. 0 means the 
%               centre is a 'peak'; pi means the centre is a 'trough',
% 


% set defaults 
if nargin < 5,         phase         = 0;
  if nargin < 4,       falloff       = 2;
    if nargin < 3,     orientation   = 0;
      if nargin < 2,   spatialFreq   = 1;
        if nargin < 1, gridSize      = 100;
        end
      end
    end
  end
end
line = linspace(-1,1,gridSize);  % create a grid of x and y coordinates, from -1 to +1
[mx,my]=meshgrid( line,line); 
z = exp(-(mx.^2+my.^2) / falloff ) .* ...
    cos(  spatialFreq  * 2 * pi  * ...
            ( sin(orientation)*mx/2 + cos(orientation)*my/2 + phase) ...
    );
  
  
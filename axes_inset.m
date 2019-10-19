function h=axes_inset(fraction, varargin)
% create a new axes as an inset to the current axes
% in the north-east corner.
%
% fraction: 0 to 1, relative size of inset axes
% 'position' : 'ne', 'nw', 'se', 'sw' == compass location of inset
%
% sgm 2015

[pos] = parsepvpairs({'position'}, {'ne'}, varargin{:});

switch pos
  case 'ne', inset=@(r,f) [r(1)+(1-f)*r(3) r(2)+(1-f)*r(4) r(3)*f r(4)*f];
  case 'nw', inset=@(r,f) [r(1) r(2)+(1-f)*r(4) r(3)*f r(4)*f];
  case 'sw', inset=@(r,f) [r(1) r(2) r(3)*f r(4)*f];
  case 'se', inset=@(r,f) [r(1)+(1-f)*r(3) r(2) r(3)*f r(4)*f];
end

h=axes( 'position', inset(get(gca,'position'), fraction) );
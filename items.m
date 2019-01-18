function varargout = items(x)
% [x1 x2 x3] = items(x)
%
% assigns each element of x to a different output of the function.
% sgm 2015
for i=1:length(x)
  varargout{i} = x(i);
end

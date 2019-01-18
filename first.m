function Y = first(X)
% helper function to get the first element of an array.
% identical to doing  X(1).
% sgm 2010
if isempty(X), Y = []; return; end
if ~iscell(X)
  Y=X(1);
else
  Y=X{1};
end
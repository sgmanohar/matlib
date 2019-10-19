% function Y=findregions(X)
% return a list of regions of true values of X
% i.e.         findregions( [1 1 1 0 0 0 1 1 1] )
% returns      [ 1, 4
%                7, 10 ]
% meaning items 1 until 4 are true, and items 7 until 10 are true.

function Y=findregions(X)
if sum(size(X)>1) > 1
  warning('findregions:vectorise', 'X is a matrix; it will be vectorised.');
end
v=diff([ 0; X(:); 0]);
Y=[];
if(isempty(v)) 
  return;
end
Y=[find(v==1) find(v==-1)];

  
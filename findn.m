function Y=findn(X)
% Y=findn(X)
% for an array X, go through each column of X and get FIND( ,1) the first 
% element that is true. Return an array of integer indices.
% e.g. if size(X) = [ 10 11 12 ], then 
%         size(Y) = [  1 11 12 ]
% note that the value is zero if no element is found.
% shortcut to get x(1). was very handy before 'find' accepted a count.
% returns [] if x is empty.
% sgm 2011

sx = size(X);
X2 = reshape(X,size(X,1),[]); % make into matrix
Y2=[];
for i=1:size(X2,2)
  w = find(X,1);
  if ~isempty(w)
    if length(w)>1 || size(Y2,1)>1 % multiple items found ?
      Y2 = nancat(2, Y2,w); 
    else % single items only
      Y2(1,i) = w; 
    end
  else
    Y2(:,i) = 0;  % nothing found
  end
end
Y = reshape(Y2,size(Y2,1),sx(2:end)); % restore to correct shape

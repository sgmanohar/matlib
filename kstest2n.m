function [h,p,ksstat] = kstest2n(X1,X2, varargin)
% 2-sample Kolmogorov Smirnov test on matrices X1, X2.
% ( see kstest2 )
% each column of X1 and X2 are samples from a different distribution.
% the results are matrices indicating pairwise comparisons of columns of X1
% and X2. P-values are NOT corrected for multiple comparisons.
% sgm 2016

% if one parameter, compare self-columns. Note this is inefficient as
% you only need off-diagonals triangle and the matrix is symmetrical...
if ~exist('X2','var'), X2=X1; end 

for i=1:size(X1,2) % for each column of X1
  for j=1:size(X2,2);% and each column of X2
    [ h(i,j), p(i,j), ksstat(i,j) ]=kstest2(X1(:,i), X2(:,j), varargin{:}); 
  end
end

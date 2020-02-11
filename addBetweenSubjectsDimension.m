
function X = addBetweenSubjectsDimension(X1,X2, X3)
%%
% X3 = addBetweenSubjectsDimension(X1,X2)
% add a between-subjects dimension at position 2, to compare two groups of
% subjects, X1 and X2.
% 
% Xi are arrays that you would normally provide to anovanTable.
% The first dimension must represent which subject.
% Each other dimesion represents a factor in the design - so Xi has the form
%  X ( SUBJECT, FACTOR1, FACTOR2, ... )
% X1 and X2 must have identical size apart from the first dimension, 
% i.e. have the same expermental design
% 
% This function inserts a new dimension 2 as a 'between subjects' factor.
% The first two dimensions of X3 then looks like
%      
%      X1_s1  nan
%      X1_s2  nan
%      ....   nan
%      nan    X2_s1
%      nan    X2_s2
%      nan     ...
% 
% so the resulting matrix has dimensions
%   size(X)  =
%     [ size(X1,1) + size(X2,1),  2,  size(X1,2) size(X1,3) ... ]
%
% and X is therefore of the form
%  
%   X ( SUBJECT, GROUP1/2, FACTORS... )
% 
% X3 is an optional third group of subjects 

S1=size(X1); S2=size(X2);
if S1(2:end)~=S2(2:end) error('dimensions of X1 and X2 dont match');
end
nd = ndims(X1);
x1=permute(X1,[1 nd+1 2:nd]); % insert singleton dimension 2
x2=permute(X2,[1 nd+1 2:nd]); 
x2=cat(2,nans(size(x2)), x2); % stretch X2 so it has two elements along dimension 2
   % which are all nan for the first element.
X =nancat(1,x1,x2);           % concatenate X1 with X2 along subject dimension,
   % such that x1 expands along dimension 2 first, giving nans for the
   % second element of dimension 2
if exist('X3', 'var')
  x3=permute(X3, [1 nd+1 2:nd]);  % add a singleton dimension 2
  s3=size(x3); s3(2)=2;
  x3=cat(2, nans(s3), x3); % stretch X3 along dimension 2 padding on the left with nans
  X = nancat(1,X,x3);
end
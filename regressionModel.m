function M = regressionModel(N, interact, nointeract)
% M = regressionModel(N, interact, nointeract)
% produce a regression model that can be used for regress, anova, x2fx etc.
% that contains all the main effects, plus all orders of interaction terms
% between the specified variables. 
%
% N = number of variables ( i.e. number of columns in the resulting model
%     matrix)
% 
% interact =  list of variables that have interactions. Each element in
%             interact must be an integer from 1 to N.
%
% e.g. N=3 and interact = [1 2 3] generates a matrix of the form
%
% [ 1 0 0 
%   0 1 0 
%   0 0 1
%   1 1 0 
%   1 0 1
%   0 1 1
%   1 1 1 ]
%   
%  so all interaction terms of the specified variables are given.
%
% alternately, use [] and specify 'nointeract' which means all variables
% will have an interaction term EXCEPT the ones in the 'nointeract' list.

if exist('nointeract','var') 
  if ~isempty(interact) error('cant use both interact and nointeract'); end
  interact = 1:N;
  interact(nointeract)=[];
end

M = eye(N);
numinteracts = 2^length(interact);
for i=1:numinteracts % for every possible combination of interaction terms
  % i is a number where each bit corresponds to a single interaction variable
  % i.e. one element of 'interact'
  c = zeros(1,N); % contrast vector
  % 'c' is the row (term) of the model that might be added
  numset = 0;
  % 'numset' keeps track of the number of bits that are set in this model term
  for j=1:N % for each variable
    if ~any(interact==j), continue; end  % if it's not an interaction variable, skip it
    jj = find(interact==j)-1; % which index in the interactions is it?
    if isempty(jj) error('whoops!') ; end
    if bitand(2^jj,i),      % does this combination i have the corresponding bit set for this variable j?
      c(j)=1; % if so, set this element of the contrast
      numset=numset+1;
    end
  end
  if numset>1,  % if we have at least a pair of interacting variables, 
    M=[M; c]; % add on a row to the model
  end
end
  


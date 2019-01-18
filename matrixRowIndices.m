function M2 = matrixRowIndices(M,I)
% matrixRowIndices(matrix, indices)
% 
% Indices is a row vector whose elements are indices  down into 
% the corresponding columns of M.
% They allow selecting by index independently, in each column of a matrix M.
%
% e.g.
% matrixRowIndices( [1  2  3  4;
%                    5  6  7  8;
%                    9 10 11 12 ], [1 3 2 2] );
% = [1 10 7 8]
%
% (i.e. choose the first item from column1, third item from column2, etc)
%
% * The indices I should be integers from 1 to size(M,1), or should be
%   boolean matrix with same size as M.
% * I can be a matrix, in which case, the result has the same number of rows
%   as I, and contains multiple selections from each row. In this case, any 
%   NaNs in I are ignored.
% * If different numbers of items are selected from each row, then the 
%   result is nan-padded.
%
% sanjay manohar 2011

if 0   % old version
M2 = M(repmat([1:size(M,1)]',1,size(M,2))==...
           repmat(I, size(M,1),1))';
         
else   % new version 2014 - also accepts booleans, or matrix of indices
  if ~islogical(I) && all(flat(isnan(unique(I)) | unique(I)==1 | unique(I)==0))
    I=abs(I)>0;  % if I is all 0s, 1s or nans, then convert to logical
  end
  sz  = size(M);                   % record original size
  szi = size(I); 
  M   = reshape(M,sz(1),[]);       % convert M into a matrix (extra dimensions flattened)
  I   = reshape(I,szi(1),[]);      % same for indices I
  M2  = [];                        % initialise result
  for i=1:size(M,2)                % go through each column of M
    tmp  = M(:,i);                 % select column of values
    tmp2 = I(:,i);                 % and column of indices
    tmp2(isnan(tmp2)) = [];        % remove nans. 
    tmp=tmp(tmp2);                 % use the indexing to extract required elements
    if isempty(tmp), tmp=nan; end; % ensure empty selections are occupied by a placeholder 
    M2=nancat(2,M2,tmp);           % and then add the selected items to result
  end
  % result should have same shape as input except for the number of rows -
  % which depends on how many items were selected by I.
  M2 = reshape(M2, [],sz(2:end));  
end
  
    

       
function a = struct2array(s, DIM)
% STRUCT2ARRAY Convert structure with doubles to an array.
% improved from the mathworks version
% sgm 2015

if isnumeric(s) && exist('DIM','var') && isstruct(DIM)
  tmp=s; s=DIM; DIM=tmp; 
end
if ~exist('DIM','var'), DIM=1; end

% Convert structure to cell
c = struct2cell(s);
if isempty(c), a=[]; return; end
if ~isnumeric(c{1}), error('not numeric'); end
if isscalar(c{1})
  a=cell2mat(c); return
end
% count dimensions of individual cells
nd=ndims(c{1});
if nd==2 && size(c{1},2)==1, nd=1; end
c=permute(c, [ [nd+1:nd+ndims(c)]  1:ndims(c) ]);

% Construct an array
a = cell2mat(c);


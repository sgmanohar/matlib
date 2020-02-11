function sout=catStruct(dim, varargin)
% syntax: 
%  1) catStruct(dim, struct1, struct2 [, struct3...])
%     concatenate two structures. Each field is concatenated along the
%     specified dimension.
%     e.g. catstruct(1, struct('a',[1 2 3],'b',[4 5 6]), 
%                       struct('a',[7 8 9],'b',[2 3 4]) )
%         = struct('a', [1 2 3; 7 8 9], 'b',[4 5 6; 2 3 4])
%  2) catStruct(dim, struct_array, fieldname)
%     concatenate the value of the field across the structure array.
%     dim is the first dimension for the field's values
%     * if 'fieldname' is '' (blank), then return all the fields combined
%       into a new structure.
%     * if the array is a cell array, use the (scalar) struct in each cell
%       element. If the cell is empty, use NaN.
% sgm 2011

if(length(varargin)<2)
  fprintf('syntax: catstruct(dim, struct1, struct2 [, struct3...])\n');
  sout=varargin{1};
  return
end

catfunc=@nancat;
s=varargin{1};
if nargin<3 || ~ischar(varargin{2}) % is third param a string? if not, 
  sout=struct();                    % do concatenation of each field
  f=fieldnames(s); 
  for(i=1:length(f)) % for each field of structure 1,
    sout.(f{i})=s.(f{i});
    for(j=2:length(varargin)) % for each remaining structure in the input
      s2=varargin{j};
      if(~isfield(s2,f{i}))
        fprintf('Warning! Structure %d has no field %s\n', j, f{i});
      else % cat the field
        q=s2.(f{i});
        sz1=size(sout.(f{i}));
        sz2=size(s2.(f{i}));
        if(length(sz1)~=length(sz2))
          fprintf('Field %s of structure %d has %d dimensions, not %d', f{i}, j, length(sz2), length(sz1));
        elseif any( find(sz1~=sz2) ~= dim )
          fprintf('Field %s of structure %d is mismatched; size is %d, not %d\n', f{i}, j, sz2, sz1);
        else
          sout.(f{i}) = catfunc( dim, sout.(f{i}), s2.(f{i}) );
        end
      end

    end
  end
else % otherwise, third param is a string. 
  f=varargin{2}; % get field name
  if isempty(f)
    if iscell(s), fn = fieldnames(s{end});  
    else          fn = fieldnames(s); 
    end
    for i=1:length(fn)
      sout.(fn{i}) = catStruct(dim,s,fn{i});
    end
    return
  end
  sout=[];
  sz=size(s);
  s=s(:);      % convert to a vector
  for i=1:length(s)
    e=s(i);
    if iscell(e), e=e{1}; end;
    if ~isstruct(e), e=struct(f,nan); end
    val = e.(f);
    if isempty(val),  val=nan; end
    sout=catfunc(dim, sout, val); % just concatenate one field
  end
  szout=size(sout); % now reshape to original size on structure dimensions
  szout(dim) = 1; % ensure the first dimension of concatenation is zero
  szout(end+1:dim+length(sz)-1) = 1; % ensure there are enough dimensions in output
  szout(dim:dim+length(sz)-1) = szout(dim:dim+length(sz)-1) .* sz; % 
  sout = reshape(sout, szout);
end
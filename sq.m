function B=sq(varargin)
% shortcyt to squeeze
% but also squeeze second dimension into first if only two dimensions
if size(varargin{1},1)==1 && length(size(varargin{1}))==2
  varargin{1}=varargin{1}';
end
B=squeeze(varargin{:});
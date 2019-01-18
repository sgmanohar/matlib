function varargout = errorBarPlotn(X,Y, varargin)
% H = PLOTN( X, Y, .... )
%
% Same as plotn, but uses errorBarPlot. it sends 3D data to the plot
% command, allowing multiple lines, and multiple plots. 
% It thus handles 5- or 6-dimensional arrays.
% Dimensions 5 are vertically arranged, and 6 is horizontal.
% Additional params recognised:
%   titlefun : a function handle, called as titlefun(i,j) where i,j are the
%              indices of the subplot (i.e. dimensions 3 and 4 of Y)
%              whose return value is used to set the title of the plot
%   eachfun  : a function handle or text string, which is evaluated for
%              each subplot. It has no input arguments - i.e. called 
%              as eachfun(), or eval('eachfun'). You can do things like 
%              set the limits of the subplot. You can use gca to get the axes. 
%   plotfun  : the plot function, which is called with a 3D value. 
%              Default is "plot", called as plot(Y), or plot(X,Y) if X given.
%              varargin is also passed to the plot command.
%              e.g. you could use "imagesc" or "contour"
%   linestyle: true if you want the 4th dimension to be line style (dots /
%              dashes etc). The function can then take 5-dimensional arrays, and
%              dimension 3 is line style, 4 and 5 are subplots. If you specify a
%              cell array here, a string from the cell array is passed to
%              plot for each line style -- e.g. {'-' ,'--', ':'} will draw
%              solid, dashed and dotted lines for dimension 3.
%   fixshape:  don't re-arrange dimensions
%   scale:     if true, then make subplot scales equal, and don't repeat
%              legends. default true.
% return: a 3- or 4-dimensional array of line handles,  similar to plot. 
% at the end of the function, ths last axes is the currently active subplot
% sgm

%%%% Parse args
if ~exist('Y','var') % only 1 param supplied
  Y=X; X=[];
elseif isstr(Y) % no X value, only arg list
  varargin=[{Y} varargin]; % Y is first element of arg list
  Y=X; X=[];
end

i=find(strcmpi(varargin, 'titlefun'));
if i
  titlefun=varargin{i+1};
  varargin(i:i+1)=[];
else titlefun=[];end

i=find(strcmpi(varargin, 'eachfun'));
if i
  eachfun=varargin{i+1};
  varargin(i:i+1)=[];
else eachfun=[];
end

i=find(strcmpi(varargin, 'linestyle'));
if i
  linestyle=varargin{i+1};
  varargin(i:i+1)=[];
  if (islogical(linestyle) && linestyle == true) || (isscalar(linestyle) && linestyle == 1)
    linestyle = {'-+',':o','--x','-.^','-v','-s','-d','-<','->'};
  elseif islogical(linestyle), linestyle=[];
  elseif isstr(linestyle) && strcmpi(linestyle,'color')
    NLS  = size(Y,3); % number of line styles
    linestyle = mat2cell(colourMap([1:NLS],NLS), ones(NLS,1),3);
  end
else linestyle=[];
end

i=find(strcmpi(varargin, 'plotfun'));
if i
  plotfun=varargin{i+1};
  if ischar(plotfun) plotfun=eval(['@' plotfun]); end
  varargin(i:i+1)=[];
else plotfun=@errorBarPlot;
end

i=find(strcmpi(varargin, 'fixshape'));
if i
  FIXSHAPE=varargin{i+1};
  varargin(i:i+1)=[];
else FIXSHAPE=false;
end

i=find(strcmpi(varargin, 'scale'));
if i
  SCALE=varargin{i+1};
  varargin(i:i+1)=[];
else SCALE=true;
end

% check data format
sz = size(Y);   % size of data
if length(sz)==3, sz(4)=1; end % in case fewer than 3 dimensions

% if the X dimension is small but there are lots of lines or subplots,
if any(sz(2:end)>9) && max(sz(2:end))>sz(1) && ~FIXSHAPE
  warning('Your data seems an odd shape. I am going to permute it to make it look nicer.\nuse "fixshape" true to prevent this in future.');
  [~,biggest]=max(sz); % find biggest dimension
  % move it to the front
  ord=1:length(sz); ord(biggest)=[]; ord=[biggest, ord];
  Y=permute(Y,ord);
  sz=size(Y);
end
if any(sz(2:end)>25) 
  warning('your data will not be comprehensible.');
end
if length(sz)<5, sz=[sz ones(1,5-length(sz))]; end % make sure 4 dimensions

%%%% do plot
if isempty(linestyle) % no linestyles; dim 2=lines, dim3/4 = subplots
  for i=1:sz(4)  % rows of subplots
    for j=1:sz(5)   % columns of subplots
      subplot(sz(4),sz(5), (i-1)*sz(5)+j);
      indx = {':',':',':',i,j}; % select 3D slice of Y.
      if ~isempty(X)
        if ndims(X)>2 % need to select slice of X also?
          x=X(indx{:});
        else x=X;
        end
        h(:,i,j)=plotfun(x,Y(indx{:}), varargin{:});
      else
        h(:,i,j)=plotfun(Y(indx{:}), varargin{:});
      end
      if ~isempty(titlefun)
        title(titlefun(i,j));
      end
      if ~isempty(eachfun)
        if isstr(eachfun), eval(eachfun);
        elseif isa(eachfun, 'function_handle'), eachfun();
        else warning('plotn:fun','unknown each-function type');
        end
      end
      if SCALE && i<sz(4) % except for ultimate row,
        set(gca,'xtick',[]); % remove x axis
      end
      if SCALE && j>1 % except for first column
        set(gca,'ytick',[]); % remove y axis
      end
    end
  end
  if SCALE, makeSubplotScalesEqual(sz(4),sz(5)); end
else % use linestyles as dimensnion 3; subplots are dim 4/5
  if length(sz)==5, sz(6)=1; end % need 5 dimensions
  if sz(4) > length(linestyle), warning('Please provide %g line styles for dimension 3',sz(3)); end
  for i=1:sz(5)  % rows of subplots
    for j=1:sz(6)   % columns of subplots
      subplot(sz(5),sz(6), (i-1)*sz(6)+j);
      for k=1:sz(4) % for each line style
        indx = {':',':',':',k,i,j}; % select 3D slice of Y
        lsi = 1+mod(k-1,length(linestyle)); % line style index
        lsargs = { linestyle{lsi} };
        if isnumeric(lsargs{1}) 
          lsargs = [{'color'}, lsargs]; 
        end
        if ~isempty(X)
          if ndims(X)>2 % need to select slice of X also?
            x=X(indx{:});
          else x=X;
          end
          h(:,i,j,k)=plotfun(x,Y(indx{:}), lsargs{:}, varargin{:});
        else
          h(:,i,j,k)=plotfun(Y(indx{:}), lsargs{:}, varargin{:});
        end
        if ~isempty(titlefun)
          title(titlefun(i,j));
        end
        if ~isempty(eachfun)
          if isstr(eachfun), eval(eachfun);
          elseif isa(eachfun, 'function_handle'), eachfun();
          else warning('plotn:fun','unknown each-function type');
          end
        end
        hold on
        % make sure colours repeat, if line plots are used
        if isequal(plotfun,@errorBarPlot), set(gca,'ColorOrderIndex',1); end
      end % next k
      hold off
      if SCALE && i<sz(5) % except for ultimate row,
        set(gca,'xtick',[]); % remove x axis
      end
      if SCALE && j>1 % except for first column
        set(gca,'ytick',[]); % remove y axis
      end
    end % next j
  end % next i
  if SCALE, makeSubplotScalesEqual(sz(5),sz(6)); end  
end % if linestyle
if nargout>1
  varargout{1}=h;
end
function varargout = plotn(X,Y, varargin)
% H = PLOTN( Y, .... )
% H = PLOTN( X, Y, .... )
%
% Same as plot, but for 3-, 4-, 5- or 6-dimensional arrays.
% plots higher dimensions as subplots on vertical then horizontal axes
% additional params recognised:
%   titlefun : a function handle, called as titlefun(i,j) where i,j are the
%              indices of the subplot (i.e. dimensions 3 and 4 of Y)
%              whose return value is used to set the title of the plot
%   eachfun  : a function handle or text string, which is evaluated for
%              each subplot. It has no input arguments - i.e. called 
%              as eachfun(), or eval('eachfun'). You can do things like 
%              set the limits of the subplot. You can use gca to get the axes. 
%   plotfun  : the plot function, which is called with a 2D value. 
%              Default is "plot", called as plot(Y), or plot(X,Y) if X given.
%              varargin is also passed to the plot command.
%              e.g. you could use "imagesc" or "contour"
%   linestyle: true if you want the 3rd dimension to be line style (dots /
%              dashes etc). The function can then take 5-dimensional arrays, and
%              dimension 3 is line style, 4 and 5 are subplots. If you specify a
%              cell array here, a string from the cell array is passed to
%              plot for each line style -- e.g. {'-' ,'--', ':'} will draw
%              solid, dashed and dotted lines for dimension 3.
%   fixshape:  don't re-arrange dimensions
%   scale:     if true, then make subplot scales equal, and don't repeat
%              legends. default true unless 'subplot' specified.
%   threed:    if true, then pass 3D array (the first 3 dimensions of Y)
%              to the plot function.  Useful for errorBarPlot, for example.
%   subplot:   false = prevent creation of subplots. Only useful for
%              when using linestyles. use @subaxis for larger panels. 
% 
% return: a 3- or 4-dimensional array of line handles,  similar to plot. 
% at the end of the function, ths last axes is the currently active subplot
% sgm

%%%% Parse args
if ~exist('Y','var') % only 1 param supplied
  Y=X; X=[];
elseif ischar(Y) % no X value, only arg list
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

FIXSHAPE = false; % prevent rearranging the dimensions for a better view?
THREED   = false; % send three dimensional slices to the plot-function?

i=find(strcmpi(varargin, 'plotfun'));
if i
  plotfun=varargin{i+1};
  if ischar(plotfun) 
    if strcmpi(plotfun,'errorbarplot')  % for errorbarplot, override the default parameters
      plotfun = @(x,varargin)errorBarPlot(x,'dostats',0,'area',1,varargin{:}); 
      FIXSHAPE = true; THREED = true; 
    else
      plotfun=eval(['@' plotfun]);
    end
  end
  varargin(i:i+1)=[];
else plotfun=@plot;
end

i=find(strcmpi(varargin, 'threed'));
if i
  THREED = varargin{i+1};  % should we pass 3D data to the plot fun?
  varargin(i:i+1)=[];
  % otherwise pass only 2D data
end

i=find(strcmpi(varargin, 'fixshape'));
if i
  FIXSHAPE=varargin{i+1};
  varargin(i:i+1)=[];
end


i=find(strcmpi(varargin, 'scale'));
if i
  SCALE=varargin{i+1};
  varargin(i:i+1)=[];
else SCALE=true;
end

i=find(strcmpi(varargin, 'subplot'));
if i
  SUBPLOT=varargin{i+1};
  varargin(i:i+1)=[];
  if ~strcmp(class(SUBPLOT),'function_handle')
    if SUBPLOT, SUBPLOT=@subplot
    else        SUBPLOT=@(i,j,k)0;
    end
  else
    % SCALE = 0; % avoid rescaling if using a subplot function?
  end
else SUBPLOT=@subplot;
end

% check data format
sz = size(Y);   % size of data
if length(sz)<4+THREED, sz=[sz ones(1,4+THREED-length(sz))]; end % make sure at least 4 dimensions

% if the X dimension is small but there are lots of lines or subplots,
if any(sz(2:end)>9) && max(sz(2:end))>sz(1) && ~FIXSHAPE
  warning('Your data seems an odd shape. I am going to permute it to make it look nicer.\nuse "fixshape" true to prevent this in future.');
  [~,biggest]=max(sz); % find biggest dimension
  % move it to the front
  ord=1:length(sz); ord(biggest)=[]; ord=[biggest, ord];
  Y=permute(Y,ord);
  sz=size(Y);
  if length(sz)<4+THREED, sz=[sz ones(1,4+THREED-length(sz))]; end % make sure at least 4 dimensions
end
if any(sz(2+THREED:end)>25) 
  warning('your data will not be comprehensible.');
end

%%%% do plot
if isempty(linestyle) % no linestyles; dim 2=lines, dim3/4 = subplots
  for i=1:sz(3+THREED)  % rows of subplots
    for j=1:sz(4+THREED)   % columns of subplots
      if strcmp(class(SUBPLOT),'function_handle') 
        SUBPLOT(sz(3+THREED),sz(4+THREED), (i-1)*sz(4+THREED)+j);
      else
        isheld=ishold; hold on; 
      end
      if ~THREED
        indx = {':',':',i,j};     % select 2D slice of Y.
      else
        indx = {':',':',':',i,j}; % select 3D slice of Y.
      end
      if ~isempty(X)
        if ndims(X)>2 % need to select slice of X also?
          x=X(indx{:});
        else x=X;
        end
        hij = plotfun(x,Y(indx{:}), varargin{:});
        if ~isempty(hij), hij=hij(1); else hij=nan; end
        h{:,i,j} = hij(1);
      else
        hij = plotfun(Y(indx{:}), varargin{:});
        if ~isempty(hij), hij=hij(1); else hij=nan; end
        h{:,i,j} = hij(1);
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
      if SCALE && i<sz(3+THREED) % except for ultimate row,
        set(gca,'xtick',[]); % remove x axis
      end
      if SCALE && j>1 % except for first column
        set(gca,'ytick',[]); % remove y axis
      end
    end
  end
  if SCALE, makeSubplotScalesEqual(sz(3+THREED),sz(4+THREED),[],[],'subplot',SUBPLOT); end
else % use linestyles as dimensnion 3; subplots are dim 4/5
  if length(sz)==4+THREED, sz(5+THREED)=1; end % need 5 dimensions
  if sz(3+THREED) > length(linestyle), warning('Please provide %g line styles for dimension 3',sz(3)); end
  for i=1:sz(4+THREED)  % rows of subplots
    for j=1:sz(5+THREED)   % columns of subplots
      if strcmp(class(SUBPLOT),'function_handle') 
        subplot(sz(4+THREED),sz(5+THREED), (i-1)*sz(5+THREED)+j);
      else isheld=ishold; hold on; end
      for k=1:sz(3+THREED) % for each line style
        try
          set(gca,'ColorOrderIndex',1);
        catch me; end % works in matlab 2016
        if ~THREED
          indx = {':',':',k,i,j};     % select 2D slice of Y
        else
          indx = {':',':',':',k,i,j}; % select 3D slice of Y
        end
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
          h{:,i,j,k}=plotfun(x,Y(indx{:}), lsargs{:}, varargin{:});
        else
          h{:,i,j,k}=plotfun(Y(indx{:}), lsargs{:}, varargin{:});
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
        if isequal(plotfun,@plot), set(gca,'ColorOrderIndex',1); end
      end % next k
      hold off
      if SCALE && i<sz(4+THREED) % except for ultimate row,
        set(gca,'xtick',[]); % remove x axis
      end
      if SCALE && j>1 % except for first column
        set(gca,'ytick',[]); % remove y axis
      end
    end % next j
  end % next i
  if SCALE, makeSubplotScalesEqual(sz(4+THREED),sz(5+THREED),[],[],'subplot',SUBPLOT); end    
end % if linestyle
if nargout>1
  varargout{1}=h;
end
if ~strcmp(class(SUBPLOT),'function_handle') && ~isheld, hold off; end

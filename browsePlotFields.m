function selectedFields = browsePlotFields(r, drawfunction, varargin)
% function browsePlotFields(r, drawfunction, ...)
%
% calls drawfunction for each field in r, but allows you to browse forward
% and backward, and search, and save figures.
%
% [pvals, anovatable] = drawfunction ( field_name ):
%   is a function that plots the given field in 'r'.
%   Usually should call plotField with the given field name.
%   It should return the p-values and anova tables if statistics are 
%   requested. This allows the browse function to search for fields which 
%   are significant for particular variables.
%   
%   if there are multiple anovas performed for this figure, then return
%   them all in a cell array of tables.
%
% Options:
%   'subplot',  [M, N] 
%      creates an array of MxN subplots, and draw all of them, pausing
%      after all have been drawn.
%
%   'filter',   string
%      The filter is the initial search pattern that is matched against the
%      fields. Only fields of r that match the filter are plotted.
% 
%      If the filter begins with '*' then the significance-line is tested
%      against the pattern - The significance line is a string shown under
%      the title which contains the factors that are significant, and
%      their p-values. For example
%           'filter', '*Age*Session' 
%      will plot only fields which have a significant interaction
%      between the factors Age and Session.
%     
%      If the filter begins with '/' then a regex match is used.
%   'saveall', 'filename_prefix'
%      saves all the figures (that match the regex), with the filenames
%   'fieldlist'  {fieldname1 fieldname2 ...}
%      just display fields with the given names


[SUBPLOT, filter, saveprefix,fieldlist] = parsepvpairs( {'subplot', 'filter', 'saveall','fieldlist'}, ...
                          {[], '', [], {}}, ...
                          varargin{:} );
SAVE_ALL = ~isempty(saveprefix);
f = fieldnames(r);
ALPHA = 0.05;
PAUSE=true;
% if this is true, then creates a new file name from the field name on each
% graph. Otherwise, use the filename from the previous save
DONT_PRESERVE_FILENAME = true;
subplotIndex = 1;
next = true;
i=1; 
if isscalar(filter), i=filter; filter=''; end % can specify a start index instead of a filter
EXIT=false;
subplotFields = containers.Map;
selectedFields = {};
fieldsPerPage = prod(SUBPLOT); if isempty(SUBPLOT), fieldsPerPage=1; end
while ~EXIT % repeat the current graph if 'next' is false
  pauseThisTime = PAUSE;
  drawThisTime = true;
  if ~isempty(fieldlist) % filter according to field list from 'fieldlist'
    if ~any(strcmp(f{i}, fieldlist))
      pauseThisTime = false;
      drawThisTime = false;
    end
  end
  if ~filterMatchField( filter,f{i} ), % filter according to current 'filter'
    pauseThisTime = false ;
    drawThisTime = false;
  else
    fprintf('MATCH ');
    fprintf('%g:%s\n',i,f{i});
  end
    
    if SUBPLOT,
      subplot(SUBPLOT(1), SUBPLOT(2), subplotIndex);
      if isnumeric(gca) % numeric axis handles (pre-2014 matlabs)
        subplotFields(num2str(gca)) = f{i};
      else              % object axis handles (post-2014)
        ax=gca;ha=DataHash(ax.Position);
        subplotFields(ha) = f{i};
      end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % call the user's draw function
    table={};
    if drawThisTime 
      [pvals table] = drawfunction(f{i});
      drawnow;
    end
    if ~isempty(table) && drawThisTime
      if ~iscell(table{1})
        disp(table(:,[1 3 6 7]));
      end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isempty(table) 
      fprintf('no statistics avaiable for field "%s"\n',f{i}); % pauseThisTime=false; 
    else
      if iscell(table{1}) % is the table an array of tables?
        pline = {};
        for k=1:length(table)
          pline = [pline titleGraph(table{k}, ALPHA, [])];
        end
      else
        pline = {titleGraph(table, ALPHA, [])};
      end
      
      if ~any(filterMatchProb(filter, pline)) % whizz by if no significant items
        pauseThisTime = false;
        drawThisTime = false; % this stops us saving the figure, if we wanted to save all without pausing, when the item doenst match.
      end
    end
    
    if SUBPLOT
      if pauseThisTime % on 'pause' items, 
        subplotIndex=subplotIndex + 1; % next subplot 
        if subplotIndex <= prod(SUBPLOT) % and if we've not reached the end of the page
          pauseThisTime = false; % carry right on to next plot
        else
          subplotIndex  = 1; % otherwise pause, and start back at top of next page
        end
      end
    end
    
    carryOn = false;
    while ~carryOn
      if(~SAVE_ALL)
        if pauseThisTime
          [x,y,key]=ginput(1);
        else
          key=[];
        end
      else key='s';
      end
      
      if any(key=='s') && drawThisTime,
        if SUBPLOT
          saveall = input('save all graphs (1/0) ?');
          if saveall==[], saveall=false; end
          if ~saveall % save selected image - enlarge to full window
            if isnumeric(gca) % numeric axis handles (pre-2014 matlabs)
              savefield = subplotFields(num2str(gca));
            else              % object axis handles (post-2014)
              ax=gca;ha=DataHash(ax.Position);
              savefield = subplotFields(ha);
            end
            clf
            drawfunction(savefield);
            drawnow
          end
        end
        if ~exist('savefield','var') || DONT_PRESERVE_FILENAME, savefield=f{i};end; % use name of last graph drawn
        fname = [saveprefix savefield '.pdf'];
        % if we're pausing after graphs, or if it's the first save, ask to
        % confirm filename.
        if pauseThisTime || ~exist('firstSaveDone','var')
          fname2 = input(sprintf('filename = %s ?',fname)); % confirm ok or change filename
          if ~isempty(fname2), fname=fname2; end
        else firstSaveDone=true; end
        % do save
        %print('-depsc', [ fname '.eps' ]);
        %saveas(gcf, fname, 'pdf');
        if  ~exist('export_fig','file') warning('cant find export_fig: please download from matworks');
        else export_fig(fname, '-pdf');
        end
        key=[]; % don't carry on - stay on same graphs
      end
      
      
      if key=='f' % change filter
        filter=input(sprintf('filter "%s" -> ', filter));
        i=i-1;
      end
      if key=='a' % change alpha
        old = ALPHA;
        ALPHA = input(sprintf('alpha %g -> ', ALPHA));
        if ~isnumeric(ALPHA) || alpha>1 || alpha<0, % validate input
          ALPHA=old; fprintf('error!  0 < alpha< 1');
        end
        i=i-1;
      end
      if key=='g'
        i = input(sprintf('go to field number (1-%g) -> ',length(f)));
      end
      if key=='h',
        fprintf('q:quit, p:previous, f:filter, c:continuous, k:keyboard, s:save image, a:alpha\n');
        i=i-1;
      end
      if key=='k',             keyboard,  end
      if key=='c',           PAUSE=false;    end
      if any(key=='q') || any(key==27),   EXIT=true;    end
      if key=='p',               i=i-2*fieldsPerPage;    end
      if i<1, i=1; end % don't go before first plot!
      % if clicked on graph, store selcted graph.
      if any(key<4) 
        if SUBPLOT
          selectedFields = [selectedFields subplotFields(num2str(gca))];
        else selectedFields = [selectedFields i]; end
      end 
      if any(key>3) || ~pauseThisTime,        i=i+1;   end
      carryOn = (~pauseThisTime) || any(key>3); % carry on with next image / screen once a key is pressed
    end
      
      if i>length(f),    EXIT=true; end
end % while not EXIT

 

function pline = titleGraph(t,ALPHA, exclude)
    % title graph based on anova table results 
    
    % get significant p values
    p=[t{2:end,7}]; sig=p<ALPHA; 
    sig(exclude)=false;
    exclude = strcmp('random',t(2:end,8));
    sig(exclude)=false;
    sig=find(sig); % list of significant terms
    ot = get(get(gca,'title'),'string'); % old title
    if ischar(ot) ot={ot}; pline=''; % make multiline
    else
      if 0 % PRESEVE OLD P-LINE?
        pline = ot{2}; % add to second line of title
      else
        pline='';
      end
    end
    for j=1:length(sig)
      pline =  [pline '; ' t{sig(j)+1,1}  ':p=' num2str(p(sig(j)))   ];
    end
    title({ot{1}, pline});
 %  text(mean(xlim),max(ylim),pline);

 
  function m=filterMatchField(filter, x)
%fprintf('matching %s == %s\n', filter, x);
    % does the field x match the filter regexp?
    if isempty(filter), m=true; return ; end
    if filter(1)=='*', m=true; return; end
    if filter(1)=='/', m=any(regexp(x,filter(2:end))); return; end
    m=any(strfind(x,filter)) ;
  
  function m=filterMatchProb(filter, x)
    % does the p-value list x match the regexp?
    if isempty(filter), m=true; return; end
    if filter(1)~='*', m=true;return; end
    if filter(2)=='/', m = cellfun( @any, regexp(x,filter(3:end)) ); return; end
    m= cellfun( @any, strfind(x,filter(2:end)) );
    %tmp = strfind(pline,filter(2:end));
    %if ~any([tmp{:}]),
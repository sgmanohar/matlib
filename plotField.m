
function [y xznames, pval, anovatable, anY, anX, h]= plotField(r, f, doAnova, plotIndividuals, factorLevels, varargin)
% [ y, xznames, pval, anovaTable, anovaY, anovaX, plothandles ] = 
%    plotField(r, f, doAnova, plotIndividuals, factorLevels, ... )
%
% Plot a field of a structure, using "CamelCase" field names.
% Structure-array r is r(subject) for each subject, with a field f.
%   i.e.  r(:).(f)
% names of fields should be as follows:
%  [ dependent_name  'For' independent_name 'And' category_name ]
% where the category is optional.
% Independent names or category names can be of the form
%  [ value1 'Or' value2 ]
% in cases where they have two possible values.
%
% an example field name could be 
%     'ReactionTimeForChooseAOrChooseB' 
%  or 'AccuracyForCueTimeAndSessionNumber'
% 
% if doAnova is true or omitted, an anova is performed on the two or three
% factors, with subject as a random factor.
% 
% plotIndividuals draws a single line for each subject, and also the
% errorbars for the group
% 
% factorLevels allows a list of names for the factor-levels to be 
% specified, which are used  for labelling the plot. They are of the form:
%    factorLevels.factorName = { 'factorLevel1', 'factorLevel2' }
% and correspond to the column/row subscripts in the data in field f. The
% factor names should be the component of the field name that is relevant.
% So if you have a field
%    r.ReactionTimeForChoiceAndSession
% then you might supply
%    factorLevels.Choice = {'ChoseA','ChoseB'};
%
% varargin: 
%  UseColourMap - useful when plotIndividuals is true
%    instead of the default "line order", create a new set of
%    colours from the COLORMAP, and gradate the colours of individual lines.
%  transform: 
%    a one-to-one function that converts a double array into a double array
%    of the same shape and size. This allows the data to be numberically 
%    transformed before plotting and statistics.
%  extraGroupings
%  extraGroupingsDimensions
%  extraVarnames
%    see anovanTable - the parameters extraGroupings and extraGroupingDimensions
%    are passed on to anovanTable directly, and the extraVarNames specify
%    the extra variable that is added by the new groupings.
%  nestExtraDimensions:
%    when extra dimensions are added, are they nested? if so, then the new
%    first of the groupingDimensions is nested within the newly created
%    grouping factor .
%  doPlot
%    if false, don't actually plot anything, just report statistics.
%  removeAllSubjectsWithNans
%    delete the rows (subjects) which contain any nans. This means that a
%    subject with even one nan will be discarded.
% 
%  continuous - a vector of integers indicating which variables are to be
%    treated as continuous
%
%  all remaining varargin{:} get sent to the errorBarPlot command, so you
%  can specify a line style for example.
% 
% return values [ y, xznames, pval ]
%   y is the full matrix of data taken from fields r.(f)
%   xznames it the parsed names of the variables, taken from 'f'
%   pval is the result of the anova, with the first element corresponding
%     to the first factor, and if the second element for the second factor
%     if present.
%   anovatable is the full table as returned by anovan (if you have asked
%     for doAnova == true). If you dont' use this return value, I will
%     print the table on the screen.

ALPHA = 0.10;
% when this is true, then dont show all the p-values: show only the
% significant ones.
TITLE_SIGNIFICANCE = true;
% if this is true, fields without a categorical axis (i.e. scalar fields
% with no 'For' in the field name) are permitted
ALLOW_BAD_FIELDS = true;

if ~exist('doAnova','var') doAnova=true; end
plotGroup = true;
if ~exist('plotIndividuals', 'var') 
  plotIndividuals=false; 
end
% plotGroup = ~plotIndividuals ; % don't plot group when plotting individuals?

removeAllSubjectsWithNans = true;
USE_COLOURMAP = true;
extraAnovanParams = {}; extraVarNames = {};
remove=[];
NO_PLOT=false;
EXTRA_GROUPINGS = false;
nestExtraDimensions = false;
transform = @(x)x;
contin = [];
SHOW_TABLE=nargout<4; % print table if not requested as output
for i=1:nargin-5
  if ~isstr(varargin{i}), continue; end
  if strcmpi(varargin{i}, 'UseColourMap')
    USE_COLOURMAP=varargin{i+1};
    remove=[remove i i+1];
  end
  if strcmpi(varargin{i},'extraGroupings') || strcmpi(varargin{i},'extraGroupingsDimension') || strcmpi(varargin{i},'nested')
    extraAnovanParams = [extraAnovanParams varargin(i:i+1)];
    remove=[remove i i+1];
    EXTRA_GROUPINGS=true; % don't allow nesting when there are extra groupings.
    extraAnovanParamStruct.(varargin{i}) = varargin{i+1};
  end
  if strcmpi(varargin{i}, 'nestExtraDimensions')
    nestExtraDimensions = varargin{i+1};
    remove=[remove i i+1];
  end
  if strcmp(varargin{i},'extraVarNames')
    extraVarNames = varargin{i+1};
    remove=[remove i i+1];
  end
  if strcmp(varargin{i}, 'doPlot')
    NO_PLOT = ~varargin{i+1};
    remove=[remove i i+1];
  end
  if strcmp(varargin{i}, 'removeAllSubjectsWithNans')
    removeAllSubjectsWithNans = varargin{i+1};
    remove = [remove i i+1];
  end
  if strcmp(varargin{i}, 'transform')
    transform = varargin{i+1};
    remove = [remove i i+1];
  end
  if strcmp(varargin{i},'continuous')
    contin = varargin{i+1};
    remove = [remove i i+1];
  end
  if strcmp(varargin{i},'showtable')
    SHOW_TABLE=varargin{i+1};
    remove = [remove i i+1];
  end
end
varargin(remove)=[];


  % parse the factor names & check for errors
  pos=strfind(f,'For'); 
  exit=false;
  if isempty(pos)  % must have at least one categorical X-variable
    warning('Badly formed field "%s" - ignoring',f); 
    if ALLOW_BAD_FIELDS 
      of=f; f=[f 'ForAll']; pos=strfind(f,'For');
      try
        [r.(f)] = deal(r.(of)); % rename the field so that it is ok
      catch mexp
        mexp
      end
    else
      exit=true;
    end
  else
    try
      tmp0=r(1).(f);
    catch me
      warning(sprintf('plotfield: No such field %s',f)); exit=true;
    end
  end
  if exit
    y=nan; xznames={}; pval=[]; anovatable={}; anX=[]; anY=[]; return; 
  end
  
  yname=f(1:(pos(1)-1)); % first part of name is Y
  xnames=f((pos(1)+3):end); % second part of name is X
  pos=strfind(xnames,'And');
  if ~isempty(pos) % is there a third part to the name?
    if length(pos)==1
      znames=xnames((pos(1)+3):end); % put it in Z
      wnames=[]; w2names=[];
    elseif length(pos)==2  % there's another category?!
      znames=xnames((pos(1)+3):(pos(2)-1)); % put them in z and w
      wnames=xnames((pos(2)+3):end);
      w2names=[];
    elseif length(pos)==3
      znames =xnames((pos(1)+3):(pos(2)-1)); % put them in z and w
      wnames =xnames((pos(2)+3):(pos(3)-1));      
      w2names=xnames((pos(3)+3):end);
    end
    xnames=xnames(1:(pos(1)-1)); % and keep only first part
  else znames=[]; wnames=[]; w2names = [];
  end
  
  nd = ndims(r(1).(f)); % how many dimensions of data?
  for(j=1:length(r)) % get each subject's data
    tmp=r(j).(f);
    if(size(tmp,1)==1 && isempty(znames)) tmp=tmp'; end;
    if exist('y','var') && (size(y,2)~=size(tmp,1) || size(y,3)~=size(tmp,2))
      error('fields dont have same size for each subject');
    end
    % Y ( SUBJ, X [, Z] )
    if(nd<3)
      y(j,:,:)=tmp; % data for this field
    elseif nd<4
      y(j,:,:,:)=tmp;
    else
      y(j,:,:,:,:)=tmp;
    end
  end
  y=transform(y);
  if all(isnan(y(:))) 
    xznames={}; pval=[]; anovatable={}; anX=[]; anY=[]; return; 
  end
  
  % do the plot
  held=ishold; h=[];
  linetypes = {'-',':','-.', '--'};
  if ~NO_PLOT
    if plotIndividuals
      if ~ismatrix(y) % plot individuals for each condition
        NI = size(y,3); % number of individuals
        if ~USE_COLOURMAP        colord=get(gca,'ColorOrder');
        else        colord=colourMap(1:NI,NI);      end
        for k=1:NI % for each z-value (each line)
          h=[h; plot(sq(y(:,:,k))', ':','Color',colord(k,:),varargin{:})];
          hold on;
        end
      else         % plot individuals - only one condition exists
        NI = size(y,1);
        if USE_COLOURMAP, colord = colourMap(1:NI,NI);
        else colord = get(gca,'ColorOrder');end
        for k=1:NI
          h=[h; plot(y(k,:), ':','Color',colord(k,:),varargin{:})];
          hold on
        end
      end
    end
    if plotGroup % plot across-subject means and errorbars
      if isempty(wnames)
        h=errorBarPlot(y,varargin{:});
      else
        for j=1:size(y,4)
          h=errorBarPlot(squeeze(y(:,:,:,j)),linetypes{j}, varargin{:});
          hold on;
        end
        hold off;
      end
    end
    if isempty(h) && ~held, hold off; end
    xlabel(deCamel(xnames));
    ylabel(deCamel(yname));
  end % if no plot

  znames2=[]; xnames2=[];
  if exist('factorLevels','var') % read from a table of factors?
    if isfield(factorLevels, xnames), xnames2 = factorLevels.(xnames); else xnames2=[]; end
    if ~isempty(znames)
      if isfield(factorLevels, znames), znames2 = factorLevels.(znames); else znames2=[]; end
      if ~isempty(wnames) && isfield(factorLevels,wnames), wnames2 = factorLevels.(wnames); else wnames2=[]; end
    else znames2=[]; wnames2=[];
    end
  end
  if isempty(xnames2) % no, we have to parse them from the variable name
    xnames2=[]; % work out the pairs of axis labels
    if(strfind(xnames,'Or')) % split at "Or"?
      orpos = [strfind(xnames,'Or') length(xnames)+1];
      xnames2={ xnames(1:(orpos(1)-1)) };
      for i=1:(length(orpos)-1)
        xnames2=[xnames2 deCamel(xnames( (orpos(i)+2):((orpos(i+1)-1)) )) ];
      end
    end
  end
  if isempty(znames2)
    znames2=[];
    if(strfind(znames,'Or'))
      orpos = [strfind(znames,'Or') length(znames)+1];
      znames2={ znames(1:(orpos(1)-1)) };
      for i=1:(length(orpos)-1)
        znames2=[znames2 deCamel(znames( (orpos(i)+2):((orpos(i+1))-1) )) ];
      end
    end
  end
  
  % title up the graph with the levels
  if ~NO_PLOT
    if ~isempty(xnames2)
      set(gca,'xtick',[1:length(xnames2)],'xticklabel',xnames2);
    end
    if ~isempty(znames2)
      hleg=legend(znames2);
      if isfield(hleg, 'title')
        set(get(hleg,'title'), 'string',deCamel(znames));
      else % new version of matlab - title is specified like this:
        try
          hleg.Title.String = deCamel(znames);
        catch
          warning('no title property on legend')
        end
      end
    elseif ~isempty(znames)
      legend(deCamel(znames));
    end
  end
  % return this cell array for the factors of anovan
  if isempty(znames), xznames={xnames}; else xznames={xnames znames}; end
  if ~isempty(wnames), xznames(end+1)={wnames}; end
  if ~isempty(w2names), xznames(end+1)={w2names}; end
  
  % statistics
  % anovaSize = 1 + (~isempty(xnames)) + (~isempty(znames)) + (~isempty(wnames));
  anovaSize = ndims(y);
  if ~EXTRA_GROUPINGS % check for nested factors - but ont if extra groupings present
    nestedAnova = zeros(anovaSize);
    for k=2:anovaSize % hunt for nested dimensions
      if all(all(all(all(sum(isnan(y),k)>=1)))), nestedAnova(1,k)=1;; end
    end
    if any(any(any(isnan(y)))) & ~any(any(nestedAnova)) % check for bad subjects
      badsubj=any(any(any(any(isnan(y),2),3),4),5);
      fprintf('%g subjects have no data! this will result in some statistics being unavailable.\n' ...
        , sum(badsubj));
      if removeAllSubjectsWithNans
        y(badsubj,:,:,:,:)=[]; % REMOVE any subject with NAN?
      else
        if 0 % PAUSE IF ERROR!
          keyboard
        else
          text(mean(xlim),max(ylim),sprintf('%g subjects missing!',sum(badsubj)));
        end
      end
    end
    if sum(sum(nestedAnova))>1
      fprintf('problem: two nestings found? ');
      keyboard
    end
    extraAnovanParams = [extraAnovanParams 'nested' nestedAnova  ];
  else % extra groupings: the extra grouping must be nested?
    if nestExtraDimensions
      nestedAnova = zeros(anovaSize+1); % one more for extra grouping
      nestedAnova( extraAnovanParamStruct.extraGroupingsDimension(1), anovaSize+1 ) = 1; % SUBJ becomes nested in the last (new) var
      extraAnovanParams = [ extraAnovanParams 'nested' nestedAnova ];
    end
  end
  extraAnovanParams = { extraAnovanParams{:}, 'continuous', contin }; 
  if isempty(znames)  % build a regression model
    varnames = {'subj', xnames extraVarNames{:}};
    model    = [1 0 ; 0 1 ] ;
  elseif isempty(wnames)
    varnames = {'subj', xnames,znames, extraVarNames{:}};
    model    = [1 0 0; 0 1 0; 0 0 1; 0 1 1 ] ;
  elseif isempty(w2names)
    varnames = {'subj', xnames,znames,wnames, extraVarNames{:}};
    model    = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; ...
                0 1 1 0; 0 0 1 1; 0 1 0 1; 0 1 1 1 ];
  else %  5 dimensions!!!
    varnames = {'subj', xnames, znames, wnames, w2names, extraVarNames{:}};
    model    = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; ...
                0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1; 0 0 1 1 0; 0 0 1 0 1; 0 0 0 1 1; ...
                0 1 1 1 0; 0 1 1 0 1; 0 1 0 1 1; 0 0 1 1 1; ...
                0 1 1 1 1];
  end
  if EXTRA_GROUPINGS
    % append a column of zeroz at the right, and add an extra row 
    % of [0 0 ... 0 1] at the end - so there are no interactions for new
    % grouping variables.
    model= [ model                    zeros(size(model,1),1);
             zeros(1,size(model,2))   1                        ];
  end
  
  if( doAnova && all(all(all(all(any(~isnan(y)))))) && ndims(y)>=anovaSize ) ...
    && ndims(y)==length(varnames), % problems in structure of ANOVA

    %if size(y,3)>5, contin=3; else contin=[]; end % continuous vars?
    %if size(y,2)>5 && isempty(xnames2), contin=[contin 2]; end   
    
    [prob, tabl, stats, terms anY anX]=anovanTable(y, 'varnames',varnames ...
        , 'random',1,'model', model ... % should I include [1 1] in the model - interaction between factor and subject
        , 'display','off' ...
        , 'ignorenan',true, extraAnovanParams{:});
    if SHOW_TABLE  % if not requesting table output, then print table 
      if 0 % ONLY SIGNIFICANT EFFECTS?
        tabl([true; false; prob(2:end)<0.15],[1 6 7])
      else % ALL EFFECTS?
        displaytable( [ tabl([1 3:end],[1 3 6 7]) ...
                        cellfun(@(x)char((x<ALPHA)*'*'), tabl([1 3:end],7),'UniformOutput',0 )]...
          ,[],[63,4,13,13,3] );
      end
    end
    pval = prob(2:end); % return the p values
    if ~any(isnan(y(:))) && any([ tabl{2:end,4} ] ) % if any singular terms in the model,
      fprintf('SINGULAR MODEL\n');
      tabl
      keyboard
    end
    anovatable = tabl;
  else
    pval = nan; anovatable = {};
  end
  if ~all(isnan(pval)) ptxt = [sprintf(' p = %g', pval)]; else ptxt=''; end
  title({ deCamel(f),  ptxt });
  if exist('tabl','var') & TITLE_SIGNIFICANCE 
    ptxt = titleGraph(tabl, ALPHA, [1]);
  end
  if ~exist('anY','var')
    anY=nan; anX=nan; 
  end

   

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

    
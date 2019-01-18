function [ R, significantFields, pvals, allPVals ] = singleSubjectStatistics(tr, gr, gn, varargin)
% [ R, significantfields, pvals, allPVals ] = singleSubjectStatistics(tr, gr, gn, varargin)
%%  WITHIN SUBJECT STATS
% tr = structure with measurement fields for each trial. Each field is a
%      vector, with one measurement for each trial.
%      Warning: All these vectors must be vertical! (i.e. along dimension 1) 
% gr = structure with grouping fields for each trial. Each grouping field
%      must be a vector with one value for each trial. The values must all
%      be integers, and must start at 1, and go up to a maximum value
%      determined by the length of the corresponding gn field.
% gn = structure with grouping names, corresponding to each field in gr.
%      Each field is a cell array of strings, one for each possible level
%      of the grouping variable. The number of levels should correspond to
%      the levels 1,2,3... of the corresponding field in gr.
% 
% This function groups trials according to each grouping field in gr, and calculates 
% the mean value of each measurable field in tr, for each condition in that grouping.
%
% e.g. if inputs are
%   tr.Rt      = rand(100,1);
%   tr.IsError = rand(100,1)>0.5;
%   gr.CueType = [ 1 2 3 2 3 1 3 2 3 1 2 ... ]' etc.
%   gr.Block   = [ 1*ones(1,50) 2*ones(1,50) ];
%   gn.CueType = {'Cue A', 'Cue B', 'Cue C'};
%   gn.Block   = {'Early','Mid','Late'};
%
% then output structure looks something like this:
%   R.meanRtForCueType      = [ 0.54, 0.45, 0.49 ]';
%   R.meanRtForBlock        = [ 0.41, 0.61 ]';
%   R.meanIsErrorForCueType = [ 0.31, 0.85, 0.40 ]';
%   R.meanIsErrorForBlock   = [ 0.55, 0.48 ]';
%
% Parameters:
%   interactions         a N x 2 cell array of field names, where each row
%                        represents a pair of grouping variables to be
%                        interacted. The resulting mean statistic looks like 
%                          R.meanRtForCueTypeAndBlock = [ a,b; c,d; e,f ]
%                        i.e. dimension 1 corresponds to levels of grouping variable 1 etc. 
%   badPairs:            a N x 2 cell array of field names, where each row
%                        represnts a measurement/grouping pair that should
%                        be omitted
%                          badPairs{:,1} are field names from tr, whereas 
%                          badPairs{:,2} are field names from gr.
%                        the field names may begin with '/' to use a
%                        regexp match on the field names to omit.
%   varianceMeasurables: calculate variance on the fields of 'tr' with these
%                        names. Cell array of field name strings.
%   medianMeasurables:   calcaulte the median for these fields in 'tr'. Cell
%                        array of field name strings
%   verbose:             give extra messages
%   plotEffects:         draw a plot for each variable that has a
%                        statistically significant result. if 2 is used,
%                        then plot All graphs irrespective of significance.
%   pause:               pause after one screenful of graphs (default=false)
%   doStats:             (default true) calculate single-subject statistics
%                        for each field - i.e. do an ANOVA across
%                        conditions for this subject -- only works for
%                        continuous variables, so the variance is
%                        meaningful!
% 
% Return values
%   R:                   The master structure of all results!
%   significantfields    a list of which fields in R had significant main effects 
%                        of the grouping factor upon the measurement, for 
%                        this single subject.
%   pvals                a list of p-values for each of the significant
%                        fields
%   allPVals             The p-value for each of the fields, significant or
%                        not, ordered in a matrix where rows correspond to
%                        measurements, and columns correspond to grouping
%                        variables.
%
% based on analyseTrioSingle
% for each measurement, analyse it according to each grouping.
% (c) sgm

[interactions, badPairs, varfields, medfields, VERBOSE, ALPHA, PLOT_EFFECTS, PAUSE, DO_STATS] = ...
  parsepvpairs({'interactions', 'badPairs', 'varianceMeasurables', 'medianMeasurables', 'verbose', 'alpha', 'plotEffects', 'pause', 'doStats'} , ...
               {{}, {}, {}, {}, true, 0.05, 1,false, true}, ...
  varargin{:});

PLOTS = [6,4]; % grid of n x n subplots
significantFields = {}; pvals=[]; allPVals=[];
f=fieldnames(tr); % trial measurement fields
g=fieldnames(gr); % grouping fields
plotindex=0;
for j=1:length(g) % check number of trials per condition
  trialsPerCondition=[];
  for k=1:length(gn.(g{j}))
    trialsPerCondition(k) = sum(flat(gr.(g{j})==k));
  end
  if VERBOSE && any(trialsPerCondition<3)
    for k=1:length(trialsPerCondition)
      fprintf('%s.\t%s \t: %g trials\n', g{j}, gn.(g{j}){k}, trialsPerCondition(k));
    end
  end
end

for i=1:length(f) % for each measurement
  y=tr.(f{i});
  if size(y,2)>size(y,1), y=y'; end;
  % if it's 2-dimensional, just leave it as it is in the results, as that
  % means it's already got a categorical aspect
  if ~isvector(y) 
    R.(['mean' f{i}]) = mean(y); 
    %continue;  % Y ( TRIAL, PREEXISTING_FACTOR  )
    if ~strfind(f{i},'For') error('no for in multidimensional field %s',f{i}); end
  end
  for j=1:length(g) % for each grouping variable
    x=gr.(g{j});
    % don't analyse if grouping field has same name as measurable!
    badPairing = strcmp(g{j},f{i});
    % check for bad field pairings, if we find a listed pairing, then omit
    % this combination.db
    for k=1:size(badPairs,1)
      if f{i}(1)=='/', fmatches = regexp(f{i},badPairs{k,1}(2:end)); else fmatches=strcmpi(f{i},badPairs{k,1}); end
      if g{j}(1)=='/', gmatches = regexp(g{j},badPairs{k,2}(2:end)); else gmatches=strcmpi(g{j},badPairs{k,2}); end
      if fmatches && gmatches, badPairing=true; end
    end
    if badPairing, continue; end % skip this pair
    
    yGrouped=[]; % store the y values for each level of the grouping
    % boolean matrix for each trial and each level of the grouping variable
    levels=1:length(gn.(g{j})); 
    levelsCheck = unique(x)'; levelsCheck(isnan(levelsCheck))=[];   % note that nans in the grouping become false, so 
    if length(levelsCheck)>length(levels), % check there are enough names 
      error('The grouping variable %s has %g levels, but only %g names specified', g{j}, length(levelsCheck), length(levels)); 
    end

    proportion=[];
    for k=1:length(levels) % for each level,
      if isvector(y)
        thisGroup = y(x==levels(k));
        if isempty(thisGroup), thisGroup=nan; end   % ensure something goes in each column
        yGrouped = nancat(2,yGrouped, thisGroup);   % put into columns according to groupings
      else % there is a preexisting factor.
        thisGroup = y(x==levels(k),:);
        if isempty(thisGroup), thisGroup=nan; end   % ensure something goes in each column
        yGrouped = nancat(3, yGrouped, thisGroup);  % YGROUPED ( TRIAL, PREEXISTING_FACTOR, NEW_FACTOR )
      end
      if i==1 % COUNTS of groupings - do this just once, not for every measurable
        proportion(k) = mean(x == levels(k));
        R.(['ProportionFor' g{j}]) = proportion;
      end
    end
    gmeans   = squeeze(nanmean(  yGrouped,1));
    gstds    = squeeze(nanstd(   yGrouped,[],1));
    gmedians = squeeze(nanmedian(yGrouped*1,1)); % 2014: *1 because, for some reason, nan/quantile doesn't like boolean input
    
    % medians??? variances???
    if isvector(y)
      combinedFname = [f{i} 'For' g{j}]; % the full field name
    else
      combinedFname = [f{i} 'And' g{j}]; % preexisting fieldnames must already include for.
    end
    R.(['mean' combinedFname]) = gmeans;
    dovar = false; for k=1:length(varfields), dovar=dovar || strcmp(varfields{k},f{i}); end; 
    domed = false; for k=1:length(medfields), domed=domed || strcmp(medfields{k},f{i}); end; 
    if dovar,    R.(['var' combinedFname])     = gstds; end
    if domed, R.(['median', combinedFname]) = gmedians; end;
    
    
    % std error  =  
    stder = gstds / sqrt(min(min(sum(~isnan(yGrouped),1)))); % std / sqrt ( number of trials in smallest group )
    stder = gstds  ./ squeeze(sqrt(sum(~isnan(yGrouped),1)));  % std / sqrt ( number of trials in that particular group )
    if PLOT_EFFECTS
      errorbar( [1:length(gmeans)],   gmeans,  stder  );
      set(gca,'xtick',[1:length(gmeans)],'xticklabel',gn.(g{j}));
    end
    % within subject stats: well, it's a 1-dimensional ANOVA across the 
    % levels of the grouping factor, unless otherwise stated. 
    t=nan;p=nan;
    if DO_STATS
      if isvector(y) && length(levels)>1 && all( sum(~isnan(yGrouped),1) > 0 )
        % as long as we have more than one level, and that all levels have at
        % least one exemplar trial,
        [p, t]=anovanTable(yGrouped, 'collapse',1,'varnames', g(j) , 'display', 'off');
      end
      if PLOT_EFFECTS
        title( {deCamel(combinedFname), ['p=' num2str(p) ]} ); ylabel(deCamel(f{i})); xlabel(deCamel(g{j}));
        drawnow
      end

      if(p<ALPHA)
        fprintf(combinedFname);
        significantFields{end+1}=combinedFname; pvals(end+1)=p;
        if (VERBOSE) t, end          % show the anova table
      end
      allPVals(i,j) = p;
    else % no stats
      p=nan;
    end
    if (PLOT_EFFECTS==1 && p<ALPHA) || PLOT_EFFECTS==2  % new plot if significant, or if plotting everything
      plotindex = plotindex + 1;
      if plotindex>prod(PLOTS), plotindex=1;
        if PAUSE,  pause; end
      end
      subplot(PLOTS(1),PLOTS(2),plotindex);
    end
  end % next grouping
  
  if isvector(y) % Interactions - only if no preexisting factors
    for j=1:size(interactions,1)
      combinedFname = [ f{i} 'For' interactions{j,1} 'And' interactions{j,2}];
      x1 = (gr.(interactions{j,1})); % the two grouping variables to interact
      x2 = (gr.(interactions{j,2}));
      % lev1 = unique(x1); lev1(isnan(lev1))=[];
      lev1 = 1:length(gn.(interactions{j,1})); % can't use unique in case there aren't all levels present in this subject
      % lev2 = unique(x2); lev2(isnan(lev2))=[];
      lev2 = 1:length(gn.(interactions{j,2}));
      proportion=[];
      for k1=1:length(lev1)   % for each combination of levels,
        for k2=1:length(lev2)
          % calculate the mean of trials that satisfy both conditions
          DIM=1;
          if size(x1,2)>size(x1,1), x1=x1'; end
          if size(x2,2)>size(x2,1), x2=x2'; end
          R.(['mean' combinedFname])(k1,k2) = nanmean( y( x1==lev1(k1) & x2==lev2(k2) ) ,DIM );
          dovar = false; for k=1:length(varfields), dovar=dovar || strcmp(varfields{k},f{i}); end;
          domed = false; for k=1:length(medfields), domed=domed || strcmp(medfields{k},f{i}); end;
          if dovar,    R.(['var' combinedFname])(k1,k2)     = nanstd(    y( x1==lev1(k1) & x2==lev2(k2) ) ,[],DIM);  end
          if domed,    R.(['median', combinedFname])(k1,k2) = nanmedian( y( x1==lev1(k1) & x2==lev2(k2) ) ,DIM);  end;
          if i==1 % COUNTS of groupings of the interaction -- only do this once, not for every measurable
            combinedCountName = [ 'ProportionFor' interactions{j,1} 'And' interactions{j,2} ];
            proportion(k1,k2) = mean(x1 == lev1(k1) & x2==lev2(k2));
            R.(combinedCountName) = proportion;
          end
        end % next k2 (next combination of levels)
      end % next k1 
    end % next interaction
  end
end
% counts


if VERBOSE, significantFields', end % display significant fields


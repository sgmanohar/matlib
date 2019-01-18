
function [betas, graphhandles, x_all,y_all] =reciprobit(latencies,flags,proportions,varargin)
% [betas, graphhandles, x_all,y_all] = reciprobit( ...
%    latencies, [flags, proportions, [ 'ParameterName', ParameterValue, [...] ]])
% plot a probit function of the reciprocals of the latencies.
% 
% 'flags' if supplied, categorises trials into several 
%    different trial types.  'flags' must be same size as latencies. 
%    Classify the latencies according to this criterion.
% 
% 'proportions' if supplied, it may determines the maximum value of the 
%         cumulative distribution. Default = 1.
%    0 :  all conditions scale to 100% (default)
%    1 :  each condition is scaled down if there is censorship 
%         (i.e. if data is outside the permitted range)
%    2 :  calculate maximum probability based upon the proportion of all 
%         trials that are used for each criterion of 'flags', and including
%         censorship.
%    a vector of size unique(flags) :
%         i.e., there should be one proportion for each level of
%         criterion, in a sorted order. This scales the proportion (y axis)
%         accordingly to the items in this vector. 
%
% ADDITIONAL PARAMS - 
%  'range' : [min max], indicates to discard saccades outside this latency 
%     range. default = [0, inf]
%  'probscale'
%  'probshift'
%     if these fine-adjustment parameters are provided, the fit
%     deviates from a true reciprocal-gaussian distribution: 
%     probscale is a fraction to decrease the width of the  gaussian by,
%     probshift is a shift in latency. 
%     default = 0
%  'plotfit' : if 0 then don't plot the best fit line. Default = 1
%  'verbose' : display information on what is being discarded / fit etc. 
%     default = 0.
%
% RETURNS: 
%  betas = [constant;  slope] of the best-fit line.
%        constant = the mean inverse latency (rate) in Hz, divided by sigma,
%          mu*sigma
%        slope is proportional to 1/standard deviation of the rate of rise, 
%          1/sigma, in time units. 
%  graphhandles = [crosses, lines], the object handles for the plotted
%        graphs, one for each criterion.
% 
% sanjay manohar

% read optional parameters
[d1 d2 range plotfit verbose legendText plotlines] = parsepvpairs( ...
    {'probscale','probshift','range', 'plotfit', 'verbose', 'legend', 'plotlines'}, ...
    {0, 0, [0,inf], 1, 0, [], 1}, varargin{:} );
if ~exist('proportions') proportions=1; end

MIN=range(1);  % ms below which to discard saccade
MAX=range(2);  % ms above which to discard saccade
ALL=1; COND=2;

if isrow(latencies)        % ensure latencies are a column 
    latencies=latencies'; 
end
if(exist('flags','var'))   % ensure flags are a column
  if isrow(flags), flags=flags';    end 
end

% If latencies is a matrix, then each column is a condition.
% So flatten the matrix into a vector and create a condition vector.
if(size(latencies,2)>1 && (~exist('flags','var') || (exist('flags','var') && prod(size(flags))<prod(size(latencies))) ))
  flags=repmat([1:size(latencies,2)], size(latencies,1),1);
  latencies=latencies(:);
  flags=flags(:);
end

% If second parameter is a matrix, then use it as the filter. 
if(exist('flags','var') && length(flags)==length(latencies) )
    filter = flags;
else              % Otherwise use all trials
    filter = ones(size(latencies));
end

%%%% Calculate Reciprobit for each condition


% Filter saccades by condition. Find unique values in the condition filter.
% Remove conditions that are 'nan'
filtervals=unique(filter); filtervals(isnan(filtervals))=[];
NF = numel(filtervals);                           % number of filter conditions
isheld = ishold(); % preserve axis state
leg    = {};       % create legend items for each filter
betas  = [];       % store slope and intercept for each filter
for fi = 1:NF                                     % for each condition filter
    fval = filtervals(fi);                        % get the condition value
    withinrange = latencies>MIN & latencies<MAX;  % use only trials in range
    selection = filter==fval & withinrange;       % select trials for this condition
    lat=sort( latencies(selection) );             % sort latencies for this condition
    if(verbose)                                   
        disp(['val=' num2str(fval) ...
          '; N=' num2str(prod(size(lat))) ...
          '; removed [' num2str(sum(filter==fval)-prod(size(lat))) '] latencies out of range']);
    end
    N = numel(lat);                              % trials in this condition
    if(N<2)
        if(verbose)
            fprintf('condition %d has insufficient trials', fval);
        end
        continue;
    end;
    maxlat=max(lat);                              % maximum latency is used to calculate 
    dt=maxlat/N; dp=1/N;                          % bin width, given number of bins N.

    mu   = mean(lat);
    sig  = std(lat);

    % divide into percentiles (already done, really, but allows
    % re-arranging the buckets)
    % can then use this for alternate 1/latency calculation
    p=[dp/2:dp:1-dp/2];
    
    % apply proportions
    if(length(proportions)==1)                   % proportion mode:
        if proportions==ALL                      % all trials:
            p=p*mean(selection);                 %  mul by proportion of all trials
        elseif proportions==COND                 % just current condition:         
            p=p*mean(withinrange(filter==fval)); %  mul by proportion in range, for current criteria
        end
    elseif length(proportions)>1                 % total supplied for each condition
        if(length(proportions)==length(filtervals))
            p=p*proportions(fi);                 % mul by individual supplied values
        else
            error('proportions must correspond to flag criteria');
        end
    end
        
    z=-ones(size(lat))./lat;                     % inverse latency

    % Treat d1 and d2 as the deviation from a normal distribution
    % probit  = sqrt(2)*erfinv(2*p-1)
    % cumdist = mu-sig*sqrt(2)*erfinv(2*p-1);
    probit    = sqrt(2)*erfinv((2-d1)*(p+d2)-(1-d1/2));

    % draw RT quantile points for each bin, for this condition
    h_points(fi) = plot(z,probit, 'Color', colourMap(fi,NF), 'Marker','.', 'LineStyle','none');
    leg{end+1} = num2str(fval); % add legend for this condition, value as specified in the category filter
    % calculate line of best fit using linear regression 
    [beta,~,resid] = regress(probit', [ones(size(probit')) z]);
    % store fit line for this condition
    betas= [betas beta];  %  betas ( intercept/slope , condition )
    hold on;
    if(plotfit)                                  % draw best fit line?
        h_lines(fi)=plot(z, beta(1)+beta(2)*z, 'Color', colourMap(fi,NF), ...
            'LineStyle', '-', 'Marker', 'none');
        leg{end+1} = [num2str(fval) ' fit'];     % legend for fit
    elseif(plotlines)                            % draw a line connecting points?
        h_lines=plot(z,probit,'Color',colourMap(fi,NF),'Marker','none','LineStyle','-');
        leg{end+1} = [num2str(fval) ' line'];    % legend for line
    else h_lines=[];                             % draw no lines?
    end
    if isempty(legendText)                       % has a legend been supplied?
      leg_crosses{fi}=num2str(fval);             % if not, use categories.
    else                                         % if so,  use it.
      leg_crosses{fi}=legendText{fi};
    end
    square_residuals(fi)=sum(resid.^2);          % calculate error term
    
    [~,ord]=sort(lat);
    x_all{fi} = z(ord);
    y_all{fi} = probit(ord);
end
if exist('x_all','var') 
  x_all=sq(nancat(x_all));
  y_all=sq(nancat(y_all));
end

if isempty(leg) % if no conditions had sufficient trials: 
    fprintf('no valid trials given');    return;
end

xlims=xlim(); xlim([xlims(1),0]);           % ensure the RT axis goes up to infinity
xvals=get(gca,'XTick'); xvals(xvals==0)=-0; % replace 0 by -0
% set up correct RT labels on the x-axis (inverse scale -1/RT)
set(gca,'XTickLabel', cellfun( @(x)sprintf('%.3g',-1./x), num2cell(xvals) , 'UniformOutput', 0));
xlabel('RT'); myerf=@(x)(erf(x)+1).*50; myerfinv=@(x)erfinv(x./50-1);
if 0,                                       % pin down the probit axis
  ylim([-3 3]);   
end   
yvals=linspace(myerf(-3),myerf(3),length(get(gca,'YTick')));
set(gca,'YTick', myerfinv(yvals) );         % and use probit proportions
set(gca,'YTickLabel',cellfun( @(x)sprintf('%.3g',x), num2cell(yvals), 'UniformOutput',0));
ylabel('cumulative probability (%)');       % add a dotted line at P = 50%
h=line([xlims(1),0],[0,0]); set(h,'LineStyle',':', 'Color',get(gca,'XColor'));
leg_crosses(h_points==0)=[]; h_points(h_points==0)=[];
legend(h_points,leg_crosses, 'Location','NorthWest');
graphhandles = [h_points h_lines];          % return the figure object handles

if(~isheld) hold off; end                   % restore the state of the axes


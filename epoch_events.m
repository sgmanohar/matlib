function Y = epoch_events(X,events, varargin)
% Y = epoch_events(X, events)
%   X is a vector time series, and events is a list of indices.
%   take the values of X( events(1):events(2) ), 
%   X( events(2):events(3) ) etc., and make a matrix where 
%   each column is a time series between a pair of events.
%   the last column is X( event(end):end ).
% 
% if X and events are matrices, then operate on each column 
% X(:,i), event(:,i), and return a 3-dimensional array.
% result = Y( timepoint, event, Xcolumn )
% 
% 'region': two-element array with [start, end]
%   if region is supplied, then each column i of y is the data in a 
%   window between event(i)+region(1) : event(i)+region(2)

i=find(strcmp(varargin,'region'));
if i>0
  region=varargin{i+1}; varargin([i i+1])=[];
else
  region = []; 
end

% for each column of X
for i=1:size(X,2)
  x=X(:,i); % get the column
  e=events(:,i); % get the corresponding list of events
  % is e a boolean vector of same length as data?
  if length(e)==length(x) && length(unique(e(~isnan(e))))==2
    % if so, treat it as flags for where events occur.
    e=find(e);
  end
  for j=1:length(e) % for each event
    % ignore nan events
    if isnan(e(j)), continue; end
    % if no 'region' parameter supplied.
    if isempty(region)
      % use the whole segment between events j and j+1
      if j<length(e) % if we're not on the last event,
        % use region from event j and j+1.
        reg = [e(j) e(j+1)];
      else  % on the last event
        % just take samples from ethe last event to the end.
        reg = [e(j) length(x)];
      end
      endnans = []; startnans = []; 
    else % if a region parameter is supplied,
      % take the area between event+region(1) and event+region(2).
      % add the event position to the region to get samples of interest
      reg = e(j)+region; 
      if reg(1)<=0 % if the left-marker is before the start,
        % pad the result with nans at the beginning
        startnans = nan(-reg(1),1);
        reg(1)=1; % and start at the beginning
      else % left-marker is within the data
        startnans=[];
      end
      % if the right-marker is after the end of the data,
      if reg(2)>length(x)
        % pad the end of the result with nans
        endnans = nan(reg(2)-length(x),1);
        % and finish at the end of the data.
        reg(2)=length(x);
      else % right-marker is within data.
        endnans = [];
      end
    end
    y{j} = [ startnans; x(reg(1):reg(2)); endnans ];
  end % next event
  Y{i} = nancat(2,y{:});
end
Y = nancat(3, Y{:});
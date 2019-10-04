function [ shiftedY, shiftedResiduals, residuals, B, dev, stats  ] = removeTimelockedRegression( Y, X, T )
%  [ shiftedY, shiftedResiduals, residuals, B, dev, stats  ] = removeTimelockedRegression( Y, X, T )
%
% Inputs: 
%   Y is an array of Y (TRIAL, TIMEPOINT_IN_TRIAL) of the value to be regressed.
%   T is the timepoints of the event of interest, as column indices into Y,
%     for each trial (so T is a column vector with as many rows as Y)
%   X is a matrix of regressors X ( TRIAL, REGRESSOR ), with a row for each 
%     trial. The regressors are used to perform linear regression at each 
%     timepoint of Y.
%
% Procedure:
% take the values of Y aligned at times T
% then perform linear regression with regressors X
%    -- by default we use GLMFIT( Y_event, X, 'normal' , 'constant','off')
% take residuals and realign them at time zero
% note: a nan value for T will make that whole trial ignored
% 
% Return values:
%   residuals: the residuals after linear regression, 
%          RESIDUALS( TRIAL, TIMEPOINT ) - i.e. it corresponds exactly to
%          the input Y.
%   
%   B:     the regrssion coefficients BETA( REGRESSOR, TIMEPOINT ) for the 
%          regression, for each time  point. 
%          it will have as many rows as there are regressors, with
%          one column for each timepoint
% 
%   dev:   these return values are from the regression, see GLMFIT
%   stats: for details
% 
%   shiftedY: the values of Y (before regression) shifted so that times are aligned on T.
%          If Y has size M x N, then shiftedY has size M x (2*N+1).
%          The points at times T are in column N.
%
%   shiftedResiduals: the regression residuals but aligned in time so that times are
%          aligned on T. Exactly the same format as for shiftedY.
%

%% Shift values according to times in T
NC = size(Y,2); % number of columns (timepoints) in Y
NT = size(Y,1); % number of rows (trials) in Y
NR = size(X,2); % number of regressors
HW = [-NC:NC] ; % time around each event to histogram
t0 = floor( T ); % ensure times are integers - for indexing 
indices = bsxfun(@plus, HW,T); % the column of Y to use for regression
badindx = (indices < 1) | (indices > NC) | isnan(indices);  % ensure indices are in-range
% INDICES will now specify the colum index of Y to use for each trial and timepoint. 
% so for example if there are 10 timepoints and T = [1; 5; 10], then 
% INDICES is [ na na na na na na na na na na  1  2  3  4  5  6  7  8  9 10 na ] 
%            [ na na na na na na  1  2  3  4  5  6  7  8  9 10 na na na na na ]  
%            [ na  1  2  3  4  5  6  7  8  9 10 na na na na na na na na na na ]
% so that timepoint 1 in trial 1, and timepoint 5 in trial 2, will be
% aligned.
indices(badindx) = 1; % temporarily substitute '1' for n/a bad indices (just to prevent problems indexing)
% create a matrix Y2 where all Y values are aligned at the corresponding
% time index from T. 
% That is, Y( 1,T(1) ) is aligned with Y( 2,T(2) ).
for tr=1:NT
  Y2(tr,:) = Y(tr,indices(tr,:));
  Y2(tr,badindx(tr,:)) = nan;
end

%% now do the regression
if nargout>1 % need to do the regression?
  Y3 = nan*ones(size(Y2));
  B  = nan*ones(NR,NC);
  for t=1:(1+NC*2) % for each timepoint
    yt = Y2(:,t); % the data for this timepoint
    if sum(~isnan(yt))>(2*NR)
      [betas dev stats] = glmfit(X,yt,'normal', 'constant','off');
      Y3(:,t) = stats.resid; % store residuals
      B(:,t)  = betas;       % and betas
      if nargout>5 % need 'stats' output? then compile values at each timepoint
        P(:,t)      = stats.p;
        dfe(:,t)    = stats.dfe;
        covb(:,:,t) = stats.covb;
        se(:,t)     = stats.se;
      end
    else
      Y3(:,t)=nan; B(:,t)=nan;
    end
  end % end for each time point
  
  if nargout>2 % do we need residuals to be realigned?
    % now realign the residuals at the original zero times
    indices = bsxfun(@minus, (NC+2):(NC*2+1), T);
    % so INDICES is now
    % [   10 11 12 13 14 15 16 17 18 19 20 ]
    % [    6  7  8  9 10 11 12 13 14 15 16 ]
    % [    1  2  3  4  5  6  7  8  9 10 11 ]
    % which gives the original indices
    % [ 0  1  2  3  4  5  6  7  8  9 ]
    % [ 0  1  2  3  4  5  6  7  8  9 ]
    % [ 0  1  2  3  4  5  6  7  8  9 ]
    badindx = (indices < 1) | (indices > (NC*2+1)) | isnan(indices);
    indices(badindx) = 1;
    for tr=1:NT
      Y4(tr,:) = Y3(tr, indices(tr,:));
      Y4(tr,badindx(tr)) = nan;
    end
    residuals        = Y4; % realigned residuals - at original timeframes
  end % don't need residuals to be realigned
  shiftedResiduals = Y3; % unre-aligned residuals (i.e. aligned at time T)
  if nargout>5 % need stats structure?
    stats.p    = P;
    stats.dfe  = dfe;
    stats.covb = covb;
    stats.se   = se;
  end
end % don't need to to regression? just shift the values
% return values
shiftedY  = Y2; 


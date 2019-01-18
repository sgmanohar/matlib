function [ t_test_result, p_vals, t_statistics, t_threshold, hplot ] ...
         = permutationOLS( Y, X, C, group, varargin )
% function [  t_test_result, p_vals, t_statistics, t_threshold, hplot ]  ...
%                  = permutationOLS( Y, X, C,  GROUP,  ['param',value]...)
%
% calculates permuted t-statistics using ordinary least squares.
% takes similar parameters to OLS.
%
% DATA       =  Y( TRIAL, SAMPLE )     
% DESIGN     =  X( TRIAL, REGRESSOR )
% CONTRASTS  =  C( CONTRAST, REGRESSOR )
% 
% results:      T_TEST_RESULT( CONTRAST, SAMPLE ) = 0 or 1
%               P_VALUE      ( CONTRAST, SAMPLE ) 
% 
% The data Y is for multiple measurements (SAMPLEs) taken on each TRIAL. 
% The statistics perform regression and calculate the threshold of t such
% that the family-wise error rate is below alpha, across the multiple
% samples.
%
% GROUP ( TRIAL ) is unique value for each permutation group (e.g.
%  subjects), such that permutations are only performed within a subject.
%  If group is omitted, all TRIALs are in the same group. 
%
% Parameters: 
%   ALPHA = false detection rate, default 0.05
%   TWO_TAILED = 0 or 1, default 1. 
%      If this is 1, then test the hypothesis "is the contrast different from zero"
%      If this is 0, then test the hypothesis "is the contrast greater than zero"
%   
% if the design is a column of ones, or is omitted, then I will 
% put -1 or +1 randomly for the permutation test, 
% to see if Y is significantly different from zero.
%
% The null distribution of (maximum-t-value over samples) is calculated,
% permuting elements that are in the same permutation group. These are then 
% compared with the t-statistics for the model OLS(Y,X,C).
%
% If P-values are requested, they are close to 0 if the contrast is significantly larger 
% than the chance distribution, and close close to 1.0 if the contrast is
% significantly lower than chance. The P-VALUE is therefore ALWAYS ONE-TAILED
% so you will need to compensate if you want a two-tailed test, i.e. use
% ALPHA/2 and 1-ALPHA/2 as p-value cutoffs.
%
% SGM 2014

if ~exist('X','var') || isempty(X)
  warning('assuming comparison of Y against zero');
  X=ones(size(Y,1),1); 
end

if ~exist('C','var') || isempty(C), C=eye(size(X,2)); end  % contrasts defaults to identity matrix
                                             % i.e. one contrast per column of X
if ~exist('group', 'var') || isempty(group),   group=ones(size(Y,1),1); end
if ~isvector(group) || size(group,1)~=size(Y,1)
  error('there should be one grouping value per row of data');
end
[alpha, TWO_TAILED, NPERMS ] = parsepvpairs( ...
  {'alpha','TWO_TAILED','NPERMS'},... % set default parameters
  {0.05   , true       , 5000   } ,varargin{:});

  
ug=unique(group);                 % grouping variables - e.g. subjects?
onesample_t =  sum(var(X))==0;    % are all rows of X the same? Then we need a 1-sample t-test.
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Estimate maximum t over all samples, 
%    from the null distribution created by within-group permutation 
t_wg_max = nan(size(C,1),NPERMS);  % this will store max and min t for each permutation
t_wg_min = nan(size(C,1),NPERMS);  
for p=1:NPERMS % for each permutation 
  pX=X;                        % initialise the permuted version of X
  for j=1:length(ug)           % then for each group:
    rows = find(group==ug(j)); % select rows for that group
    if ~onesample_t            % standard t test   -- H0: "gradient == 0"
      pX( rows, : ) = pX( rows(randperm(length(rows))), :); % re-arrange just those rows
    else                       % one-sample t-test -- H0: "value == 0"
      pX( rows,: ) = ((rand(size(rows))<0.5)*2-1) .* pX( rows,: ); % if so, then choose -1 or +1 randomly
    end
  end
  % now calculate within-group permutation's statistic:
  % T ( CONTRAST, SAMPLE, PERMUTATION )
  [~,varb,t_wg] = ols(Y,pX,C); % T_WG ( CONTRAST, SAMPLE )
  
  % calculate the distribution of max( t ) across all samples
  % so that we can correct for the family-wise error rate of finding a
  % single false positive across the whole range of samples.
  t_wg_max(:,p) = max(t_wg,[],2);   % maximum T across samples: T_WG_MAX ( CONTRAST, PERMUTATION )
  if TWO_TAILED
    t_wg_min(:,p) = min(t_wg,[],2); % minimum T across samples (for 2-tails)
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. calculate the T corresponding to a FDR of alpha, for each contrast..
%    t_statistics ( CONTRAST, SAMPLE )
[beta,~,t_statistics] = ols(Y,X,C); % do least squares with the true regressors

if TWO_TAILED
  upper_t_threshold = quantile(t_wg_max',1-alpha/2); % work out cutoff values of t corresponding 
  lower_t_threshold = quantile(t_wg_min',  alpha/2); % to the required family-wise error rate 
  t_test_result = bsxfun(@gt, t_statistics, upper_t_threshold') ...
                - bsxfun(@lt, t_statistics, lower_t_threshold');       % put +1 or -1 if contrast is above or below zero
  t_threshold  = cat(1,upper_t_threshold, lower_t_threshold); % T_THRESHOLD ( UPPER/LOWER, REGRESSOR )            
  % For how many permutations is the maximum (across samples) permuted t-statistic 
  % greater than the true t-statistic (for each sample)? If it is rarely
  % above the real t, then p is small, because t is bigger than expected by
  % chance.
  p_vals_upper = mean( bsxfun(@gt, t_wg_max , permute(t_statistics, [1,3,2])) , 2  ); % mean across permutations,
  % How often is the permuted minimum t-statistic across samples lower than 
  % the true t-statistic of each sample? If it is rarely below the real t,
  % then p is small, because t is lower than expected by chance.
  p_vals_lower = mean( bsxfun(@lt, t_wg_min , permute(t_statistics, [1,3,2])) , 2  ); % for each contrast and sample
  p_vals       = min( p_vals_upper, p_vals_lower );  % give whichever p-value is lower. Correct for 2-tails.
  % the t_test_result will tell us which direction it it significant in.
  % The returned p-value is close to 0 if either tail is significant, but
  % must be tested with a threshold of ALPHA/2. 
else % ONE TAILED
  t_threshold   = quantile(t_wg_max',1-alpha);         % t cutoff values for given FWER 
  % is the t-statistic bigger than alpha-percent of the permuted max-t
  % statistics?
  t_test_result = bsxfun(@gt, t_statistics, t_threshold');        % is the statistic bigger than the threshold?
  % what proportion of the permuted maximum-t statistics are above the true
  % t-statistic for each sample?
  p_vals        = mean( bsxfun(@gt, t_wg_max , permute(t_statistics, [1,3,2])) , 2  ); 
end
% P_VALS (CONTRAST, SAMPLE)
p_vals = permute(p_vals,[1,3,2]);                    
if iscolumn(p_vals) p_vals=p_vals'; end

if nargout>4 % Reqeusted hplot? 
  % PLOT the t test results. This shows the significant points where
  % delta is nonzero, below the graph as a bar.
  % Note: this might upset adding any more subplots to the figure!
  mainaxes = gca;
  rect=get(mainaxes,'position'); % select region in bottom 5% of axis
  axes('position',[rect(1) rect(2) rect(3),rect(4)*0.05]);
  hplot=imagep(p_vals); axis off; colorbar off;
  axes(mainaxes); % revert to main axes
end

  
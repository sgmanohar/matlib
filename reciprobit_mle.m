function [beta t0 nll]=reciprobit_mle(RT, COND, varargin)
% reciprobit_mle(RT, COND, varargin)
% fit reciprobit by mle
% RT = reaction times - vector.
% COND = conditions. if not specified, conditions can be given as columns
%        of RT, with RT as a matrix. NaN is ignored.
% return: 
%  BETA ( SIGMA/MU, CONDITION )
%  T0   = delay
%  NLL  = negative log likelihood of model


[range]=parsepvpairs({'range'}, {[0 inf] }, varargin{:});

% confine RT to range
if length(range)~=2, error('range must have 2 values, upper and lower bound'); end
RT(RT<range(1) | RT>range(2))=nan;

% convert condition vector into a matrix of RTs, if specified
if exist('COND','var') && ~isempty(COND)
  if ~isvector(RT) warning('column structure of RT will be ignored'); end
  conds = unique(COND); conds(isnan(conds))=[]; NC=length(conds); RT2=[];
  for i=1:NC; f=COND==conds(i); RT2=nancat(2,RT2, RT(f)); end; RT=RT2;
end % now RT is a matrix, each column is one COND.

% probability of RT=t given params t0 (delay), s=variability in rate,
% m=mean rate of rise.
% fn_prob (t, t0, s, m) = normpdf( (1./(t-t0) - m)./s );
% expanding singletons, to allow parameters to be vector.
fn_prob = @(t, t0, s, m)   normpdf( 1./(t-abs(t0)), repmat(m,size(t,1),1) , repmat(s,size(t,1),1) );
% parameters are given as sigma and mu for each condition, and a single
% value of t0 across all conditions.
% p ( t0, s1, s2, ..., sn, m1, m2, ..., sn )
fn_prob_trial = @(p) fn_prob(RT, p(1), p(2:1+NC), p(2+NC:1+2*NC) );
% sum over trials and conditions to get the full log likelihood
fn_loglike    = @(p) nansum(nansum(log(eps+fn_prob_trial(p))));

sig0 = 0.004; mu0 = 0.005; t00 = nanmean(RT)/10;
p0 = [0.03, repmat(sig0,1,NC), repmat(mu0,1,NC)];
% [p1, nll]=fminsearch( @(p)-fn_loglike(p), p0 , optimset('MaxFunEvals',3000));
[p1, nll]=fminsearchs( @(p)-fn_loglike(p), ...
  [0 0 0 0 0 0 0;0.1 0.1 0.1 0.1 0.1 0.1 0.1] , 10, optimset('MaxFunEvals',3000));
beta = reshape(p1(2:end), NC, 2)'; % beta ( SIG_MU, COND )
t0   = abs(p1(1));
if any(beta(:)>1e3 | beta(:)<-1e3)
  keyboard
end

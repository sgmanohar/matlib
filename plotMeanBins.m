function r=plotMeanBins(x,criterion, subjectcriterion, normalise,varargin)
% r=plotMeanBins(x,criterion, subjectcriterion, normalise,varargin)
%
% plot means of values in x, binned according to criterion and subject, 
% plotting between subject error bars.
% if 'normalise'=1, then plot the main effect of criterion - i.e.
% subject means subtracted from values before taking criterion means and sd.
% 'normalise'=2 uses multiplicative normalisation i.e. ratio to subject
% mean.
%
% x, criterion, and subjectcriterion are arrays with same size.
% the subsequent arguments are passed to 'plot'.

ucr=unique(criterion);
usc=unique(subjectcriterion);
if(~exist('normalise','var'))  normalise=0; end
for(i=1:length(ucr))
    for(j=1:length(usc))
        c=(criterion==ucr(i)) & (subjectcriterion==usc(j));
        m(i,j)=nanmean(x(c));
        s(i,j)=nanstd(x(c));
    end
end
% mean across subjects
ms=nanmean(m,1);          % mean of each subject across criteria
mc=nanmean(m,2);          % mean of each criterion across subjects
if(normalise==1)
  mn=bsxfun(@minus, m,ms); %subtract subject means
elseif normalise==2
  mn=bsxfun(@rdivide,m,ms); % normalised data
else
  mn=m;
end

mm=nanmean(mn,2);            % mean of each criterion
sm=nanstd(mn,0,2);           % std across subjects for each criterion

if(normalise==1)
  mm=bsxfun(@plus,mm,mc);  % add on criterion mean across subjects
elseif(normalise==2)
  mm=bsxfun(@times,mm,mc);  % rescale means
  sm=bsxfun(@times,sm,mc);  % rescale SDs
end

r.criteria=ucr';
r.mean=mm';
r.std=sm';
r.data=m;
errorbar(r.criteria, r.mean, r.std/sqrt(length(usc)),varargin{:});
xlabel('condition'); ylabel(['mean +/- SE across ' num2str(length(usc)) ' subjects']);

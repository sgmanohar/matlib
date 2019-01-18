function h=meanReciprobit( RT, CRITERIA, NQ, varargin)
% meanReciprobit ( RT{SUBJECT}(TRIAL) , CRITERIA{SUBJECT}(TRIAL), NUM_QUANTILES );
% returns a handle to the plot
% extra params after NQ  get sent to PLOT.

if ~exist('CRITERIA','var') || isempty(CRITERIA)
  CRITERIA=cellfun(@(x)true(size(x)),RT,'uniformoutput',0); 
end
if(~exist('NQ','var') || isempty(NQ)) NQ=50; end
for(i=1:length(RT)) % subject
  uf=unique(CRITERIA{i}); uf(isnan(uf))=[];
  for(j=1:length(uf)) % criterion
    f=CRITERIA{i}==uf(j);
    rt=RT{i}(f); %rt(rt<50 | rt>800)=nan;
    lastbin=-inf;
    for(k=1:NQ) % quantile
      q=quantile(rt,(k-0.5)/NQ);
      mrt(i,j,k)=q;
      lastbin=q;
    end
  end
end
probit=erfinv( ([1:NQ]-0.5)/NQ * 2 - 1 );
holding = ishold();
for(j=1:length(uf))
  leg{j}=sprintf('%d',uf(j));
  % if no params specified, use colormap colours
  if isempty(varargin) ,  % 1 = override color
    extraparams = {'color',colourMap(j,length(uf))};
  else extraparams = {}; end
  h=plot(-1./sq(nanmean(mrt(:,j,:),1)), probit, varargin{:}, extraparams{:});
  hold on;
end
if ~holding, hold off; end
legend(leg);
set(gca,'yticklabel',(erf(get(gca,'ytick'))+1)/2,...
        'xticklabel',-1./get(gca,'xtick') );

   
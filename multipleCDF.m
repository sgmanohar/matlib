function [P,support] = multipleCDF(valueCells, legnd, varargin)
% multipleCDF(valueCells, legnd)
% takes a series of cells containing vectors of values,
% and plots a CDF for each cell.
c=colormap;
nc=size(colormap,1);
SMOOTHING= 5;
for i=1:length(valueCells)
    if(length(find(~isnan(valueCells{i})))<2) 
        warning(['criterion ' num2str(i) ' contains no data']);
        continue;
    end;
    [f0,x0]=ecdf(valueCells{i});
    if any(strcmpi(varargin,'smooth')), 
      f0=smooth(f0,SMOOTHING,'moving'); x0=smooth(x0,SMOOTHING,'moving'); 
    end
    plot(x0,f0,'Color', c(floor((i-1)*nc/length(valueCells))+1,:));
    hold on;
    leg{i}=num2str(i);
    X0{i}=x0;
    F0{i}=f0;
end;
if(exist('legnd') && ~isempty(legend))
    legend(legnd);
elseif exist('leg')
    legend(leg);
end;
hold off;

if nargout>0
  allX0 = vertcat(X0{:});
  support = linspace( min( allX0 ), max(allX0), 100 );
  for i=1:length(X0)
    x0=X0{i}; f0=F0{i};
    if length(unique(x0)) ~= length(x0) % any nonunique values?
      [~,q]=unique(x0); % find unique ones only
      x0=x0(q); f0=f0(q);
    end
    bad = isnan(x0) | isnan(f0);
    x0(bad)=[]; f0(bad)=[]; 
    P(:,i) = interp1(x0,f0, support,'nearest','extrap' );
  end
end

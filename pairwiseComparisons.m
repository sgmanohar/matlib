function [varargout] = pairwiseComparisons(X, varargin)
% function [p, t] = pairwiseComparisons(X)
% X ( SUBJECTS, CONDITIONS )
% 
% runs a paired t-test between thc conditions represented by columns of X.
% 
% Returns CONDITIONS x CONDITIONS matrixes of p-values and t-values.
% Note that only the upper triangle will contain values - to avoid
% duplicate comparisons.
% 
% if the output of the statistical function is a vector, then return a
% 3-dimensional array.
% 
%  statfn - a function of two vectors. The outputs of this function will be
%           returned in matrix form.
%  draw_p - draw horizontal lines for significant pairs. 
%           assumes the X-axis has values 1,2,3... corresponding to the
%           columns of X ( e.g. as drawn by "errorBarPlot(X)" )
%           If a custom statfn is provided, the value should specify which 
%           number output of the statfn corresponds to p-values.
%           e.g. for @ttest, p is the second output, so use 'draw_p',2.
%  alpha  - p-value threshold for drawing lines
%  text   - the format string for p-values (see sprintf). default = '%0.2g'
% e.g. 
%  [ h, p, ci, stat ] = pairwiseComparisons( rand(100,3), 'statfn', @ttest2 )  
% 
%  h  = 
%       NaN    0    0
%       NaN  NaN    0
%       Nan  Nan  Nan
%  p  = [3 x 3 matrix]
%  ci = [3 x 3 x 2 matrix]
%  stat = 
%       []   [1x1 struct]  [1x1 struct]
%       []             []  [1x1 struct]
%       []             []            []
%      
[statfn DRAW_P , ALPHA, TEXT]=parsepvpairs( ...
  { 'statfn',    'DRAW_P' , 'ALPHA' ,'TEXT'    },...
  { @test   ,    []       , 0.05    ,'p=%0.2g' } ...
  ,varargin{:}  );

if isstr(statfn), switch statfn
    case 'ttest2',    statfn=@test2; % these functions return the p-value 
    case 'ttest',     statfn=@test;  % as the first output! 
    case 'signrank',  statfn=@signrank;
    case 'ranksum',   statfn=@ranksum;
    otherwise      error('unknown function %s',statfn);
end;end

% create a blank output cell array
N = size(X,2); % number of variables to compare
out = cell( 1, abs(nargout(statfn)) );
for k=1:length(out), varargout{k}=nan(N); end 

for i=1:N
  for j=i+1:N % for each pair of columns of X, 
    [out{:}] = statfn(X(:,i),X(:,j)); % run the statistics
    for k=1:length(out) % for each output of the statistics, 
      if isvector(out{k}) && isnumeric(out{k}) % put it in a matrix if it's a number
        % varargout{k} = nanassign( varargout{k}, [i,j,nan], out{k} ) ;
        if size(varargout{k},3)~= numel(out{k}),  
          varargout{k}=cat(3,varargout{k}, nan(size(varargout{k},1), size(varargout{k},2), numel(out{k})-1) ); 
        end% increase dimensions if needed
        varargout{k}(i,j,:) = out{k};
      else % otherwise put it in a cell.
        if isnumeric(varargout{k}), varargout{k}=cell(N); end
        varargout{k}{i,j} = out{k};
      end
    end
  end
end


if DRAW_P
  p=varargout{DRAW_P};
  yl=ylim; vpos = yl(2)-0.05*(yl(2)-yl(1)); dv = (yl(2)-yl(1))*0.025;
  nsig=0;
  for i=1:N
    for j=i+1:N
      if p(i,j)<ALPHA
        nsig=nsig+1;
        y = vpos + (nsig-1)*dv;
        line([i j], [y y]);
        line([i i], [y y-0.5*dv]);
        line([j j], [y y-0.5*dv]);
        if TEXT
          txt=sprintf('p=%0.2g',p(i,j));
          text( (i+j)/2, vpos+nsig*dv ,  txt ,'fontunits','normalized','fontsize',0.05);
        end
      end
    end
  end
end

function [p,t] = test(x,y)
[~,p,~,t] = ttest(x,y);
t=t.tstat;

function [p,t] = test2(x,y)
[~,p,~,t] = ttest2(x,y);


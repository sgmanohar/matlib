function multiplot(Y,X,L,R,C, CATEGORICALX)
% multiplot(Y,X,L,R,C, CATEGORICALX)
%
% Visualise up to 5 dimensional data!
% plot the values in Y, or a set of X/Y pairs, as a function of 3 other
% categorical variables.
% 
% each matrix must have same dimensions (or have size [1,1], in which case
% there will be only one item on that axis).
% Y = the data values themselves, which will be averaged together and
%     plotted as the y-coordinate
% X = the x coordinate for each corresponding data point: should be
%     categorical variable if CATEGORICALX is true.
% L = the line that each point belongs to - i.e. which series
% R = which row of subplots each point belongs to
% C = which column of subplots each point belongs to
% 
% L, R or C may be omitted (or []).
%
% Basically this function is a 5-dimensional plot command.
% Let's say you have one measurement Y, and it varies as a function of
% 4 categorical variables, X1, X2, X3 X4. You could call
%    multiplot( Y, X1, X2, X3, X4, true )
% and the different values of X1 appear on the X-axis, 
%         different values of X2 appear as different coloured plot-lines
%         different values of X3 appear as rows of subplots
%         different values of X4 appear as columns of subplots
% 
% 
% sgm

if (exist('CATEGORICALX')~=1) 
    CATEGORICALX=1;
end;


colours='rbgcymk';
colour=@(x)colours(mod(floor(x),length(colours))+1);

if ~exist('L','var') || isempty(L) || isscalar('L'),  L=ones(size(Y));end;
if ~exist('R','var') || isempty(R) || isscalar('R'),  R=ones(size(Y));end;
if ~exist('C','var') || isempty(C) || isscalar('C'),  C=ones(size(Y));end;

uc=unique(C)'; nc=length(uc);
ur=unique(R)'; nr=length(ur);
ul=unique(L)'; nl=length(ul);
maxysc=-inf; minysc=inf;
for c=1:length(uc)
    for r=1:length(ur)
        subplot(nc,nr,c+(r-1)*nc);
        hold off;
        for l=1:length(ul)
            x=X(C==c & R==r & L==l);
            y=Y(C==c & R==r & L==l);
            if(CATEGORICALX)
                xc=unique(X)';
                yq=bsxfun(@eq, xc, x);
                yrep=repmat(y,1,length(xc));
                ydat=yrep.*yq;
                ymu=nanmean(ydat,1);
                yse=nanstd(ydat,0,1)/sqrt(length(xc));
                plot(1:length(xc),ymu, [colour(l) 'x-'], ...
                     1:length(xc),ymu-yse, [colour(1) '.:'], ...
                     1:length(xc),ymu+yse, [colour(1) '.:']);
            else
                scatter(x,y, 3, [colour(l) '.']);
            end
            hold on;
        end
        hold off;
        if(CATEGORICALX)
            xc=unique(X)';
            lab=num2str(xc');
            for(i=1:size(lab,1))label{i}=lab(i,:);end;
            set(gca, 'XTick', 1:length(xc));
            set(gca, 'XTickLabel', label);
        end;
        mm=get(gca, 'YLim');
        if(maxysc<mm(2)) maxysc=mm(2);end;
        if(minysc>mm(1)) minysc=mm(1);end;
    end;
end;
for c=1:length(uc)
    for r=1:length(ur)
        subplot(nc, nr, c+(r-1)*nc);
        set(gca, 'YLim', [minysc, maxysc]);
    end
end;
                


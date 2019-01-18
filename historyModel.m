function [b bint r rint stats names M]=historyModel(Y,X, varargin)
% [b bint r rint stats names regressorMatrix] = historyModel(Y,X, [,'useY'] [,regression options])
% this works like a multiple regression [ regress( Y,X ) ], but tries to 
% explain Y in terms of X and the 'history of X'.
% 
% Y is a series of event outcomes, in time order, and factors that influenced
% previous events can also effect future events.
% X is a matrix of factors (regressors) affecting the corresponding value
% of Y (as well as future values of Y)
%
% useY: flag to allow the previous values of Y to themselves influence 
%       future values of Y. Default=0.
% noBlank: flag to prevent adding in a blank ('column of ones') regressor
%       if one is not supplied
% names: name each regressor; the following parameter should be a cell
%       array of strings the same length as X's regressor dimension
% history: the length of history to use. Each factor's history is used as
%       a regressor for H future items. default=1
% x:    add interaction term. each interaction is a triplet containing a
%       factor number (i.e. column number of X), its history position, and
%       a second factor number and its history position. i.e. 
%       [factor1, history1, factor2, history2]
% plot: draw a scatterplot of the result of regression
% includeonly: includes only the flagged items in the actual regression
%       calculation. The excluded items are left in place for calculation
%       of history effects  - i.e. they appear in correct position under 
%       the (trial N-1) effects etc. This is useful for regressing over
%       a random subset of trials, then predicting over the remainder.
% quiet: don't print any output. default 0 (prints output)
%
% the output 'names' gives the list of the names of the regressors
% the output 'regressorMatrix' gives the final regressor (including X) with
%     the appropriate history values. If an 'includeonly' subset is
%     specified, the regressorMatrix still includes all excluded items.
%     However it is shorter than the original Y, as the first H items are
%     missed off from Y.
% SGM 2012

% check for options
useY=0; noBlank=0; remove=[]; plotit=0;
include=ones(size(Y)); % items to include
ix=[]; % interactions
quiet=0;
for i=1:length(varargin)
    if strcmpi(varargin{i},'useY') 
        useY=1;
        remove=[remove i];
    end
    if strcmpi(varargin{i},'noBlank')
        noBlank=1;
        remove=[remove i];
    end
    if strcmpi(varargin{i},'names')
        names=varargin{i+1};
        remove=[remove i i+1];
        i=i+1;
    end
    if strcmpi(varargin{i},'history')
        H=varargin{i+1};
        remove=[remove i i+1];
        i=i+1;
    end
    if(strcmpi(varargin{i},'x'))
        nix=varargin{i+1};
        if(size(nix,2)~=4) error('interaction parameter x should be 4 items across'); end;
        ix=[ix;nix];
        remove=[remove i i+1];
        i=i+1;
    end
    if(strcmpi(varargin{i},'plot'))
        plotit=1;
        remove=[remove i];
    end;
    if(strcmpi(varargin{i},'includeonly'))
        include=varargin{i+1};
        if(size(include)~=size(Y)) error('include filter must be same size as Y'); end;
        remove=[remove i i+1];
        i=i+1;
    end
    if(strcmpi(varargin{i},'quiet'))
        quiet=1;
        remove=[remove i];
    end;
end
varargin(remove)=[];
if(~exist('names','var')) names={};
end;
if ~isempty(ix) && (any(any(ix(:,[2,4])))>H)
    error('interaction items must not have history further back than the regressors');
end;

% Align/rotate data
if(size(X,1)~=size(Y,1))
    % flip data if not correctly aligned
    if size(X,2)==size(Y,2) 
        X=X';  Y=Y';
    else
        error('X and Y must have the same length');
    end;
end

for(i=1:size(X,2))
    if(i>length(names)) names{i}=sprintf('Factor %d', i);
    end;
end;


if(~exist('H','var'))
    H=1; % history length
end;
%%%% make historical regressors
M=X;
for(i=1:H)
    M=[M(2:end,:), X(1:end-i,:)];
    for(j=1:size(X,2))
        names=[names, sprintf('%s(N-%d)',names{j},i)];
    end;
    if(useY)
        M=[M, Y(1:(end-i), :)];
        names=[names, sprintf('Y(N-%d)',i)];
    end
end;


for (i=1:size(ix,1))
    h1=ix(i,2); h2=ix(i,4);
    c1=X(H-h1+1:end-h1, ix(i,1));
    c2=X(H-h2+1:end-h2, ix(i,3));
    M=[M,  c1.*c2 ];
    names=[names, sprintf('%s(N-%d)x%s(N-%d)', names{ix(i,1)}, ix(i,2), names{ix(i,3)}, ix(i,4)) ];
end;
    
%%%% check for a blank regressor
if(~any(all(X==ones(size(X)), 1)))
    if(~noBlank)
        M=[M, ones(size(M,1),1)];
        names=[names, 'Const'];
    end;
end;

%%%% check covariance
cm=abs(triu(nancov(M),1));
[i,j]=find(cm>0.5);
for(k=1:length(i))
    if(~quiet)
        %fprintf('Factors %d and %d have covariance %d\n', i(k),j(k),cm(i(k),j(k)));
    end
end;


%%%% perform regression
Y=Y((1+H):end);         % only predict the later trials (first ones don't have antecedent history
Yf=Y(include((1+H):end) ~= 0);     % limit to included items
Mf=M(include((1+H):end) ~= 0,:);       % limit regressor to included lines
[b bint r rint stats]=regress(Yf, Mf, varargin{:});

[i,i]=sort(-abs(b));
if(length(i)>50) i=i(1:50);end;
if(~quiet) % output the largest valued regressors and their beta coefficients
    for(j=1:length(i)) 
        if prod(bint(j,:))>0 sig='*';
        else                 sig=' ';
        end;
        fprintf('%s\t%g\t%s\n', names{i(j)}, b(i(j)), sig); 
    end;
end;

%%%% plot?
if(plotit)
    scatter(Mf*b,Yf,'.');
    xlabel('X*beta'); ylabel('Y');
    held=ishold();
    hold on; plot([0,max(Y)], [0,max(Y)]);
    if(~held) hold off; end;
end;


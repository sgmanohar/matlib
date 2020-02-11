function [ m ] = nanmeanAllArgs( varargin )
% nanmeanMeanAllArgs - take nanmean of each argument
m=[];
for(i=1:length(varargin))
    m=[m, nanmean(varargin{i})];
end
end


function [ m ] = nanstdAllArgs( varargin )
% nanmeanMeanAllArgs - take nanmean of each argument
% horizontally concatenates the values for nanStd of each argument.
% sgm
m=[];
for(i=1:length(varargin))
    m=[m, nanstd(varargin{i})];
end
end


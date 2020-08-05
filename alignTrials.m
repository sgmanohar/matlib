function [d,s] = alignTrials( d, s )
% scan through a experiment-generated result structure d
% (created using transpIndex(result.data) 
% and in parallel scan through a Eyelink-generated saccade structure s
% (created using readAllEdfAsc() or readEdfAsc() )
% and for each trial in the experiment, locate the appropriate saccade.
% it corrects for duff trials and practice trials.

j=0; % counter for saccades
ok_js =  []; % keep a list of the ok saccades
for i=1:length(d.R)
  problem = 1; 
  startj = j;
  while problem && j <= length(s)% keep incrementing the saccade pointer j until we are ok
    j=j+1; % next saccades
    bts = textscan( s(j).BT_m, '%d T %d') ;
    b_s = bts{1};
    t_s = bts{2}; % get trial and block from s
    if b_s == d.block(i) && t_s == d.trialIndex(i) 
      problem = false; % we foudn the right saccades
    end
  end % next saccades
  if problem % we didn't find anything matching
    error('no trial found %g',i);
  end
  ok_js = [ok_js j]; % add the current saccades index to our list
end

s=s(ok_js);


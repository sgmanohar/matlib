
  function nstrings=deCamel(strings)
  %    new_strings = deCamel( old_strings )
  % Remove "camel case", putting spaces, to make a more friendly, readable
  % name. For example,
  %  'meanRewardSensitivity' --> 'mean Reward Sensitivity'
  %
  % Can take a cell array of strings, or a single string.
  % Useful for converting field names or file names into readable labels.
  % sgm 2011
  
  single=~iscell(strings);
  if single 
    strings={strings};
  end
  for si=1:length(strings)
    s=strings{si};
    prevcap=false; ns=[]; 
    for i=1:length(s)
      if s(i)==upper(s(i)) & isletter(s(i))
        if prevcap % another capital? just add it
          ns=[ns s(i)];
        else % no prev cap? 
            ns=[ns ' ' s(i)];  % prev word was >1 charput a space first
        end
      else % lower case (or nonletter): just add it
        ns=[ns s(i)]; 
      end
    end
    ns(ns=='_')=' '; % replace underscores with spaces too.
    nstrings{si}=ns;
  end
  if single
    nstrings=nstrings{1};
  end
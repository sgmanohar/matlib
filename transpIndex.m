function str=transpIndex(struc, ignorefieldnames)
% struct = transpIndex ( struc , [ ignorefieldnames ] )
%
% shift the dimensions of a struct array into the dimensions of its fields
% works with maximum two-dimensional struct-arrays.
% 
%
%%% NOTE:
% 1-dimensional fields only (at present!)
% e.g.
%  x = 6x3 struct array with fields
%      fa: [1x1 double]
%      fb: [5x1 double]
%
%  transpIndex(x) =
%      1x1 struct array
%      fa: [6x3 double]
%      fb: [6x3x5 double]
%
% ignorefieldnames should be a cell of field names to leave out of the
% transposed array. 
%
% sanjay manohar

if prod(size(struc))<4
    warning('you are trying to transpose indices of a struct-array smaller than 4 items')'
end;
if length(struc)==0
  str = struct; % empty struct
  str(:)=[];    % really very empty
  return;
end
fnm=fieldnames(struc);
str=struct();
ignorefield=false(length(fnm),1);
if exist('ignorefieldnames','var')
    for i = 1:length(ignorefieldnames)
        ignorefield = ignorefield | strcmp(fnm, ignorefieldnames{i});
    end
end
for i=1:size(struc,1)
    for j=1:size(struc,2)
        for k=1:size(fnm)
            if ignorefield(k) continue;end;
            dt=struc(i,j).(fnm{k});
            if islogical(dt) dt=1*dt; end; % convert logicals to numbers (allowing nans)
            if (isnumeric(dt) | islogical(dt) | ischar(dt)) & prod(size(dt))>1
                tm=prod(size(dt));
                if isfield(str, fnm{k})
                    ttm=size(str.(fnm{k}),3);
                    if tm>ttm % if trial longer, add a new slice of nans in all other trials
                        for m=ttm+1:tm
                            str.(fnm{k})(:,:,m)=nan(size(str.(fnm{k}),1), size(str.(fnm{k}),2));
                        end;
                    end;
                end;
                for m = 1:tm % put data
                    str=setfield(str, {1}, fnm{k}, {i,j,m}, dt(m));
                end;
                ttm=size(str.(fnm{k}),3);
                if ttm>tm    % if shorter, pad this trial with nans
                    for m=tm+1:ttm
                        str=setfield(str,{1},fnm{k},{i,j,m}, NaN);
                    end
                end;
            else
                %fprintf('set str.%s(%d,%d)=',fnm{k},i,j);
                % disp(dt);
                if prod(size(dt))==0 % empty
                  if( (~isfield(str, fnm{k}) || ~isstruct(str.(fnm{k})) ) )
                    try
                      str=setfield(str, {1} ,fnm{k}, {i,j}, NaN);
                    catch e
                    end
                  else % it's a field already, but it's a structrure, can't put nan in it
                    %str(i).(fnm{k})(i,j) = struct(); % leave blank...
                  end
                else
                    if iscell(dt) 
                        warning(['ignoring ' fnm{k} ': contains cells.']);
                        ignorefield(k)=1;
                    elseif isstruct(dt)
                        warning(['ignoring ' fnm{k} ': contains structs.'])
                        ignorefield(k)=1;
                    else
                        str2=setfield(str, {1} ,fnm{k}, {i,j,1}, dt);
                        str=str2;
                    end;
                end;
            end;
        end;
        fprintf('.'); % each trial
    end;
end;

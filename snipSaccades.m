function [r, info]=snipSaccades(s, start, finish, varargin)
% [A, info]=snipSaccades(saccadeStructure, start, finish [,'param', value...])
% 
% snip out saccade paths for each trial.
%
% With no extra parameters, simply returns the nan-padded array of saccades
% as complex numbers. 
% Also handles pupil size, rotating/reflecting, baseline/centre subtraction.
%
% input: 
%   saccadeStructure s: is a struct-array created by the program readEDFASC
%             and contains the field pos(t,1:4) = [times, x,y,pupil]
%             and fields with the names in 'start' and 'finish', containing 
%             times of the start and finish events.
%   start:    A string specifying the field name in the structure to use as
%             start times for cutting.  If 'start' is a number rather than 
%             a string, begin at this fixed datapoint index. If it is a
%             vector (the same size as s), begin at this index into the 
%             respective trial
%   finish:   as for start, except absolute numbers and vectors indicate 
%             the length to include after the start point. i.e., 
%             If 'finish' is a number rather than a string, take exactly this 
%             number of data points after 'start'.
%
% parameters:
%  'saccadeonly', N : only retain the eye position data for the actual
%             time of Nth saccade.
%  'rotate', X : specify a vector the same size as s, with the angle X to rotate
%             each saccade, in radians. Default 0.
%  'flip', B : a vector the same size as s, of boolean values indicating
%             whether or not to flip each trial. If no flip angle is
%             provided, the reflection is left-to-right. Default 0.
%  'flipangle' X : a vector the same size as s, of angles in radians (or an
%             [n x 2] matrix of vectors) representing the axis along which
%             to reflect each saccade.
%  'centre' [x,y] : the screen centre in pixels, to be subtracted from all x and
%             y coordinates of pos. Can have one row, or one row per trial.
%             Can be a complex number instead of [x,y]
%  'stretch', L : interpolate the snipped saccades so that each trace is 
%             interpolated to have length L, and so they all line up to
%             the given length L. i.e., there is no padding of nans at the
%             end of the saccade - it is stretched in time to fit.
%  'baselinesub', X: subtract baseline values. the parameter can be: 
%             X = single number N: then the mean baseline is taken over N 
%                     samples after event1
%             X = string 'event3': the baseline is taken at the point event3
%  'pupil' :  analyse pupil data instead of eye position. In this case, the
%             values are all reals.
%  'meancriteria', X : take the mean values of the eye data, across trials.
%             The parameter X specifies which group/bin each trial falls
%             into. The return value has one row for each unique value of
%             the criteria, in sorted order. info contains the criteria in 
%             order.
%  'smoothing', N : average using a boxcar function over N samples 
%  'plot' :   plot the values on a graph. Any unrecognised input parameters 
%             are sent to the plot command.
%  'clean', T :  number of samples for cleanup. It is the minimum length
%             of valid samples needed to retain a segment of data. Default
%             50 samples (50 msec).
%  'excludenan': exclude trials where any datapoint is nan, from the mean
%             and from the plot. (Does not exclude them from the returned
%             traces -- only from means!).
%  'minsize' : minimum saccade size to register, in pixels. Use 1 degree 
%             (~30 pixels) to remove "microsaccades". 
%  'saccadeonly' N : get rid of all data except the actual position data
%             for the Nth saccade(s)
%  'verbose' : display extra information about what is being calculated
%  'browse' : for each trial, display a diagram of the trace of the saccade. 
%             1: skip through fast, 2: pause for key, 3: enter debugger
%  'traces' : if 0, then don't return an array of traces. [] is returned as
%             the first return value. Speeds up calculation of saccades,
%             but is not compatible with plot / browse / meancriteria
%             option.
%  'interpolate', N: if a sequence of fewer than N nans occurs, after
%             cleaning, then use linear interpolation to fill in the gap.
%  'extrapolate', N: if the epoch begins or ends with a sequence of fewer 
%             than N nans, then replace them with the first/last valid value. 
%  'hampel', N: hampel ("median deviation from the median") filtering
%             removes large deviations from median caused by sample noise.
%             The filter uses a window of the specified half-width. 0 means
%             no filtering. Default 0.
%  'speedfilter', 0/1: remove (flatten) velocities which are greater than a
%             threshold. This is useful for removing sudden vacillations in
%             pupil size. Note that this will remove saccade artifacts, but
%             also remove saccades if used on the eye position! default 0.
% returns: 
%   a matrix with one row for each trial, containing the snipped eye data
%

%%%% parse parameters
[rotate, flip, centre, flipangle, stretch, baselinesub,...
 pupil, meancriteria, CLEAN_MS, saccadeonly, minsize,...
 verbose, browse, excludenan, doplot, RETURN_TRACES, ...
 smoothing, interpolate, extrapolate, HAMPEL, speedfilter]...
 = parsepvpairs( ...
    {'rotate', 'flip', 'centre','flipangle', 'stretch', 'baselinesub',...
     'pupil', 'meancriteria', 'clean', 'saccadeonly', 'minsize', ...
     'verbose','browse','excludenan','plot', 'traces', ...
     'smoothing', 'interpolate', 'extrapolate','hampel','speedfilter'}, ...
    { zeros(size(s)), zeros(size(s))>0, [0 0], 1j, 0, 0, ...
      0, [], 50, 0, 0, ...
      0, 0, 0, 0, 1, ...
      3, inf , inf , 0, 0 }, varargin{:} ...
    );

if size(flipangle,2)==2 flipangle=angle(flipangle*[1,1j]); end; % convert to column of complexs
if length(flipangle)==1 flipangle=flipangle*ones(size(s)); end; % duplicate it if single value given
if(stretch==1) stretch=100;end; % default as 100
if(length(meancriteria)==1) meancriteria=meancriteria*ones(size(s)); end;
if ~RETURN_TRACES & (browse | doplot | meancriteria )
    error('no traces is not compatible with plotting / mean options');
end

if(stretch>0) warning('off', 'MATLAB:interp1:NaNinY'); end;

% other global parameters that you can fiddle with:
MS_PER_SAMPLE      = mean(diff(s(1).pos(:,1))); % sampling rate
VSMOOTH            = 8;     % milliseconds for smoothing of velocity.
ALLOW_TRANSGRESS   = false; % include saccades that cross the start- or end-times specified?
SMOOTH_IN_CLEANING = 3 ;    % apply a minor smoothing while cleaning, so that points around bad points are excluded.
MAX_INTERP_DIST    = 1000;  % maximum difference between two samples, to permit interpolation
smoothMethod       = 'moving'; % use 'boxcar' for spktools version of 'smooth'
%% process saccades
%%%% results for each saccade:
r=[];
info.sRT    = []; % onset time of saccade
info.sEndpt = []; % physical location of endpoint
info.sBlink = []; % is there a blink in this saccade? 1/0
info.sAmpl  = []; % amplitude = abs(sVec)
info.sBendT = []; % time of maximal deviation from saccade line
info.sBendA = []; % angle of point with maximal deviation from saccade line (but now distance of this point from the saccade line)
info.sVec   = []; % endpoint relative to startpoint
info.sCurvA = []; % unsigned maximal curvature (radial acceleration)
info.sCurvS = []; % signed maximal curvature 
info.sDepA  = []; % angle of departure (where it leaves the 30-pixel circle)
info.sSpd   = []; % maximum speed of eye during saccade
info.sDur   = []; % duration of saccade

prev_warn_state = warning('off','nancat:emptyskipped');
for(i=1:length(s))  % for each trial
          scStart=nan; scEnd=nan; % start with blank values for this trial
        scBegpt=nan; scBlink=nan; scBendT=nan; scBendA=nan; scCurvA=nan;
        scCurvS=nan; scAmpl=nan; scCurvT=nan; scEndPt=nan; scBendP=nan;
        scBendR=nan; scDepA=nan; scSpd=nan;  

  if ~isempty(s(i).pos) % are there any samples?
    %%%% Epoching
    %%%% (Snip from start to end)
    if(isnumeric(start))  % start time can either be:
        if(length(start)~=1)  t1=start(i);       % vector of start positions
        else                  t1=start;          % fixed start position
        end
    else                      t1=s(i).(start);   % event name
    end
    if(isnumeric(finish))
        if(length(finish)~=1) t2=finish(i)+start;% vector of finish positions
        else                  t2=finish+t1;      % fixed finish position relative to start
        end
    else                      t2=s(i).(finish);  % event name
    end
    if ~isempty(t1) && ~isempty(t2) % ensure there are times
    
    filter = s(i).pos(:,1) > t1 & s(i).pos(:,1) < t2; % select rows of 'pos' by their timestamp.
    
    pt=s(i).pos(filter,1);                     % times of points in region of interest
    if(pupil==0)                               % eye position is columns [2,3]
        p =s(i).pos(filter,[2:3]) * [1;1j];    % get whole trace for this saccade
        if size(centre,1)==length(s)           % allow per-trial centres
          ctr = centre(i,:);
        else ctr = centre;
        end
        if size(ctr,2)==2, ctr=ctr*[1;1j]; end % convert to complex number, if not already
        p=p-ctr;                               % centre trace on zero
        p=p*exp(rotate(i)*1j);                 % rotate trace
        if(flip(i))                            % reflect along axis angled at 'reflectangle' radians
            p=abs(p) .* exp( 1j * (2*flipangle(i)-angle(p)) );
        end
    else
        p=s(i).pos(filter,4);                  % pupil trace is column 4
    end                                        
    original_p = p; % keep track of raw data. p will get post-processed!
    % NOW: p is either pupil size or position (complex), depending on 
    % whether 'pupil' = true.  "pt" are the corresponding sample times.
    if sum(filter)<100
      warning('Segments less than 100 ms in trial %g',i);
      %continue; 
    end
    
    %%%% calculate speed
    if ~exist('smooth','file')                 % smoothing function - does it exist?
      warning('cannot find smooth function');
      % create 'smooth' function in the current workspace... but don't
      % allow it to think smooth is a variable!
      eval('smooth=@(x,y,z)nanconv(x,ones(y,1)/y,''same'');'); 
    end
    v=diff(p);               % dp/dt - derivative of position or pupil
    spd=abs(v);              % speed = absolute value of derivative
    if(length(v)>VSMOOTH)    % as long as there are more than NSMOOTH timepoints, 
      vs=smooth(v,VSMOOTH,smoothMethod); % then apply some smoothing, using an NSMOOTH-millisecond filter
    else vs=[];              % if the trace is very short, then don't calculate velocity
    end
    if ~pupil % this section applies only to saccades!
      if ~isempty(s(i).saccade) % check there some saccades were recognised by eyelink
        %%%% Measure the saccades
        as=diff(vs);            % grab points where vel > 3 px/ms
        % flag to indicate if saccade is during the period of interest
        if ALLOW_TRANSGRESS     % allow saccades that transgress start and endpoints?
          saccin=s(i).saccade(:,2)>t1 & s(i).saccade(:,1)<t2; 
        else                    % don't allow transgressing start and endpoints
          saccin  =s(i).saccade(:,1)>t1 & s(i).saccade(:,2)<t2; 
        end
        scStartT=s(i).saccade(saccin,1); % times of start (for each trial)
        scEndT  =s(i).saccade(saccin,2); % and end of saccades
        
        
        if(length(scStartT)>0 && length(scEndT)>0)  % if saccades aren't all empty:
          for(j=1:length(scStartT))                 % for each saccade recorded,
            scStart(j,:) = find(pt>=scStartT(j),1); % sample index of start
          end % now, sometimes the saccade doesn't terminate correctly.
          for(j=1:length(scEndT)) % for each corresponding end-time,
            if scEndT(j)> pt(end) % if the parsed saccade end time is after the the actual trace
              scEnd(j,:) = length(pt); % make the saccade end at the end of the trace
            else                  % otherwise the saccade ended normally, 
              scEnd(j,:) = find(pt>=scEndT(j),1); % so find the real end point.
            end
          end
          
          % now calculate the saccade vector: from start-point to end-point
          % direct calculation: 
          %   scVec  =( s(i).saccade(saccin,[6,7])-s(i).saccade(saccin,[4,5]) ) * [1;1j];
          % But here, use transformed coords, from 'p':
          scVec  = p(scEnd)-p(scStart);
          if(size(scVec,2)==1) scVec=scVec'; end; % ensure it's a row-vector
          % remove small saccades
          remove = abs(scVec) < minsize;
          scStartT(remove)=[]; scEndT(remove)=[]; scVec(remove)=[];
          scStart( remove)=[]; scEnd( remove)=[];
          % now just check in case something went wrong. There should be 
          % equal numbers of saccades in the 'start' and 'vec'. 
          % So this should never happen....
          if(length(scVec)~=length(scStartT)) disp(scStartT, scVec);keyboard
          end
          
          % now we have located the saccades in the samples, 
          % we can extract the trajectories
          for j = 1:length(scStart)  % for each saccade
            scPath    = p(scStart(j):scEnd(j)) - p(scStart(j)); % the path of the sacc, relative to point of 
            scBegpt(j)= p(scStart(j)); % beginning point 
            scEndpt   = scPath(end);   % vector of this saccade
            scEndPt(j)= p(scEnd(j));   % endpoint of each saccade
            scBlink(j)= any(isnan(scPath)); % will catch nans in x or y
            scAmpl(j) = abs(scEndpt);  % length of the vector from start to end
            % calculate deviation from a straight line, at each p
            % i.e. let e be the unit saccade direction (scEndpt "hat")
            % then deviation = position minus e*(position "dot" e).
            scDev     = scPath - scEndpt.*(real(scPath).*real(scEndpt) ...
              + imag(scPath).*imag(scEndpt))./(abs(scEndpt).^2);
            % where and when is the maximum deviation from straight line?
            maxDevPt  = find(abs(scDev)==max(abs(scDev)),1); % first point of max deviation
            scDevMax  = scDev(maxDevPt); % maximum vector of deviation (it's perpendicular to the trajectory).
            if ~isempty(maxDevPt) % && ~any(isnan(p(maxDevPt-1:maxDevPt+1))) % make sure surrounded by decent datapoints! 
              scBendT(j)= maxDevPt;
              scBendA(j)= mod(pi+angle(scPath(scBendT(j))) - angle(scEndpt),2*pi)-pi;
              scBendR(j)= abs(scDev(maxDevPt)) ... actual deviation distance, 
                * (2*(imag(scDevMax * scPath(maxDevPt)'.' )>0)-1); % make it the same sign as imag( dev . traj* ), which is dev dot perp(traj)
              scBendP(j)= p(scBendT(j));
            else % we could not find a point of maximum deviation
              scBendT(j)=nan; scBendA(j)=nan; scBendP(j)=nan;
            end
            if length(scPath)>5 % if the saccade path is at least 5 samples long,
              % calculate saccade velocity, over time, using the mean of 
              % 3 consecutive velocity measures (i.e. using 4 samples)
              scVel     = diff(smoothn(scPath,3,smoothMethod)) / MS_PER_SAMPLE;
              % calculate this on *unit* velocity if you want d/ds (curvature) rather than d/dt
              scAcc     = [0;diff(smooth(scVel,3,smoothMethod))] / MS_PER_SAMPLE;
              % projection of acceleration onto unit vector perpendicular to velocity
              % ddp - (ddp . dp) dp / |dp|^2
              %     = Radial component of acceleration
              scCurv    = scAcc - (real(scAcc).*real(scVel) + imag(scAcc).*imag(scVel)) .* scVel./(abs(scVel).^2);
              sccurvt   = find(scCurv==max(scCurv),1);
              if ~isempty(sccurvt),    % saccade has a max-curvature time:  
                scCurvT(j) = sccurvt;  % the time at which the curvature (radial acceleration) is greatest
                scCurvs   = abs(scCurv) .* (real(scCurv).*imag(scEndpt) - real(scEndpt).*imag(scCurv));
                scCurvA(j)= scCurvs(scCurvT(j));
                scCurvS(j)= scCurv(scCurvT(j));
              else % problem with this saccade's curvature?
                scCurvT(j) = nan;   scCurvA(j) = nan; scCurvS(j) = nan;
              end
              % calculate the sign term as the cross product of the endpoint and the nontangential acceleration
              % find "point of departure", when saccade actually gets moving! 
              departT   = find( abs(smooth(scPath,3,smoothMethod)) > 50,1 ); % arbitrary 30 pixel criterion
              % calculate max speed: peak velocity. 
              % disallow movement of > 100px in 1 ms.
              scSpd(j)  = max(abs(scVel)); if scSpd(j)>100, 
                badpos = find(scSpd(j));
                if length(badpos)==1 && badpos==length(scSpd)
                  scSpd(j) = max(abs(scVel(1:end-1)));
                  if scSpd(j)>100, scSpd(j)=nan; end
                else
                  scSpd(j)=nan;
                end
              end
              
              % depature angle: the direction the saccade goes in, at the
              % start of its trajectory (i.e. when it leaves a 50 pixel
              % radius of the start point)
              if ~isempty(departT)   scDepA(j) = angle( scPath(departT) ); 
              else scDepA(j)=nan;  % only valid for larger saccades.
              end
            end
          end
        end
      end % if empty(saccades) . 

      % return saccade segment only? cut out all of data except Nth saccade
      if(saccadeonly && any(~isnan(scStart)))
        if length(scStart)>=saccadeonly % ensure we have enough saccades to get the Nth one
          p=p(scStart(saccadeonly):scEnd(saccadeonly)); % extract positions
          v=v(scStart(saccadeonly):scEnd(saccadeonly)-1); % and velocities
          spd=spd(scStart(saccadeonly):scEnd(saccadeonly)-1); % and abs(vel)
        else
          p=[]; v=[], spd=[]; % we don't have enough saccades on this trial
        end
      end


    end % end if ~pupil. --- Next section applies both to saccades and pupil
    
    %%%% Cleaning
    % Removes regions around blinks.
    % blinks are identified by NaN values in the samples, but extend to
    % include regions around the NaNs during which the samples are not
    % stable. 
    if ~isempty(v) && ~saccadeonly % 2014 added no cleaning if saccadeonly is requested
      % Around blinks, remove samples where speed greater than BLINKSPD
      % was 1.5
      
      nanC=nan+nan*1j; % complex nan: both x and y are nan
      if(pupil) 
        % What is the maximum absolute value of the sample that is permitted?
        % pupil sizes > 20,000 and eye positions > 2000 pixels are excluded.
        ABSMAX   = 20000;
        % speed above which data is removed around nans. Normally 1.5 or
        % 2.5
        BLINKSPD = 3.5;    
      else
        ABSMAX   = 2000;
        BLINKSPD = 0.5;    
      end;
      % find regions in which the speed is too high too.
      bad=find(   (spd>40 & [abs(diff(v)); 0]>30) ...
                | (abs(p(1:end-1))>ABSMAX) | (abs(p(2:end))>ABSMAX)  );
      p(bad)      = nanC; % turn these bad samples into NaN
      p(isnan(p)) = nanC; % note that if real(p) is nan, it doesn't guarantee imag(p)==nan.
      % now scan through the whole trace, looking for NaNs
      lastnan=1; lastwasfast=0;
      removedsegments(i)=0;
      if SMOOTH_IN_CLEANING   % this makes samples surrounding NaN
        % also into NaN - expanding the "bad" region a little.
        spdsm= smooth(spd,SMOOTH_IN_CLEANING,smoothMethod);  
      else             % the speed calculated here is used to determine
        spdsm= spd;    % "stability" of the trace around the NaN.
      end              % This helps get a good velocity when there is sample noise.
      % remove short segments of data surrounded by nans.
      % scan forward and backward from the nan, in range CLEAN_MS, and
      % remove points where the speed abs(spdsm) > BLINKSPD
      for k=1:length(p)         % for each sample, 
        if isnan(p(k))          % if it is NaN,
          if ((k-lastnan)<CLEAN_MS) && ((k-lastnan)>1)  % if there was a recent NaN,
            p(lastnan:k)=nanC;   % then erase the trace from the previous NaN until now.
            removedsegments(i)=removedsegments(i)+1;
          end
          lastnan=k; lastwasfast=0;
          % remove fast segments before nans. Go back sample-by-sample,
          % looking at the speed. If the speed remains higher than the
          % BLINKSPD, then erase the samples.
          m=1;while( k-m>0 && ~isnan(p(k-m)) && ...
              ( spdsm(k-m)>BLINKSPD || (k-m-1>0 && spdsm(k-m-1)>BLINKSPD) || ...
              (k-m-2>0 && spdsm(k-m-2)>BLINKSPD) ) ...
              );
            p(k-m)=nanC; %  Keep going until two consecutive non-fast samples have been confirmed.
            m=m+1;
          end
        else %  not nan: remove fast movements after nans
          % if the previous sample was fast or was NaN, and this or the next
          % sample is fast, 
          if ( (k-lastnan==1) || lastwasfast) && k<length(p) && ...
             ( spdsm(k)>BLINKSPD || (k+1<length(p) && spdsm(k+1)>BLINKSPD) || ...
              (k+2<length(p) && spdsm(k+2)>BLINKSPD) );
            p(k)=nanC;     % then erase this sample
            lastwasfast=1;
          else
            lastwasfast=0; % otherwise end the "fast" sequence
          end
        end
      end;
      
      %%%% Pupil filtering
      % these pupil-specific filters help remove partial-blinks and saccade
      % artifacts in the data.
      if pupil    % if we are examining pupil data,
        if speedfilter
          % Speed filter simply sets high speeds to zero. This means that
          % areas of steep change become perfectly flat. Uses BLINKSPD as a
          % "speed limit".
          v=diff(p);             % calculate speed at each instant;
          vs=smooth(v,10,smoothMethod); % averaged over about 20 ms
          v(abs(vs)>BLINKSPD)=0; % set speed to zero if it's too big;
          v(isnan(v))=0;         % nan samples have zero speed;
          firstValidP = p(find(~isnan(p),1)); % value of first valid sample 
          if isempty(firstValidP) firstValidP=nan; end % perhaps the trial is all nans?
          isnanP = isnan(p);     % keep track of invalid samples
          % Now re-integrate the new speeds to get new position
          p = firstValidP + [0;cumsum(v)]; 
          p(isnanP)=nan;         % restore invalid samples
        end
        if HAMPEL % hampel filtering? 'hampel' is the half-width of the filter, in samples
          p=hampel([1:length(p)]',p,HAMPEL,3);
        end
      end % if pupil

      %%%% Interpolation
      if interpolate
        startnans=1; lastwasnan=0; % keep track of sequences of nans
        for k=1:length(p)               % for each sample
          if isnan(p(k)) && ~lastwasnan % is this the beginning of a sequence of nans?
            startnans=k;                % if so, store the current index
          elseif lastwasnan && ~isnan(p(k)) % is this the end sequence of nans?
            if   (k-startnans <= interpolate) ... % is the sequence of NaNs short enough?
                && (startnans>1) ...    % no extrapolation!
                && (k<length(p)-2) ...  % no extrapolation if only 2 points at end
                && (~isnan(p(k+1))) ... % ensure not extrapolating to single point 
                && (~isnan(p(k+2))) ... 
                && abs(p(startnans-1)-p(k))<MAX_INTERP_DIST % and ensure distance isn't too big
              % then do the interpolation for this segment of NaNs
              p(startnans:(k-1)) = interp1([startnans-1,k], [p(startnans-1), p(k)], [startnans:(k-1)]);
            end
          end
          lastwasnan = isnan(p(k));
        end
      end
      %%%% Extrapolation
      % This just looks at the beginning and end of the epoch, and if there
      % are nans all the way (e.g. if the epoch started during a blink),
      % then assume the value of the first non-nan sample in this period.
      % Same thing applies to missing data at end of epoch.
      % Extrapolation is necessary becase initial value is needed for 
      % baseline subtraction! 
      if extrapolate 
        firstValid = find(~isnan(p),1);      % Find index of First non-nan sample.
        if firstValid>1                      % If first sample wasn't valid,
          p(1:firstValid-1) = p(firstValid); % then overwrite initial nan-segment
        end                                  % with the first valid value.
        lastValid = find(~isnan(p),1,'last');% Find Last non-nan sample.
        if lastValid<length(p)               % If the last sample wasn't valid,
          p(lastValid+1:end) = p(lastValid); % then overwrite final nan-segment
        end                                  % with the last valid value.
      end
    else % not doing cleaning - v=[] or saccadeonly=true
      removedsegments=[];
    end
    hasnan(i)=any(isnan(p)); % does the trial's trajectory contain any NaNs?
    
    %%%% Trimming
    % remove saccades that have only a nan segment after them.
    if ~pupil && ~isempty(s(i).saccade) &&  ~isempty(scStart) && ~saccadeonly% we have nonempty saccades?
      for j=1:length(scStart)                                 % go through each saccade
        if isnan(scEnd(j)) continue; end                      % some might be duds
        % remove saccades where pretty much all the points are nan after
        nextgood = scEnd(j) + find(~isnan(p(scEnd(j):end)),1); % next good point after saccade
        later    = min(length(p),nextgood+500);     % Previously: find a point 500 ms later (no longer used)
        if isempty(nextgood) ||  ... % if no samples after this saccade, or
          (nextgood<=length(p) && mean(isnan(p(scEnd(j):nextgood)))>0.95) || ... % if 95% of future samples are NaN
          scEnd(j)==length(p), % if the saccade ends at or after the end of the trial,
        
          scStartT(j)=nan; scVec(j)=nan; scEndpt(j)=nan; scBlink(j)=nan;
          scAmpl(j)=nan; scBendT(j)=nan; scBendR(j)=nan; scCurvA(j)=nan;
          scCurvS(j)=nan; scDepA(j)=nan; scSpd(j)=nan; scEndT = nan;
          %info.sVec(end,j)=nan; info.sCurvA(end,j)=nan; info.sCurvS(end,j)=nan; info.sEndpt(end,j)=nan;
          %info.sBendT(end,j)=nan; info.sBendA(end,j)=nan;  info.sAmpl(end,j)=nan; info.sDepA(end,j)=nan;
        end % then reject this saccade
      end
    end
    
    %%%% Smoothing
    if smoothing & length(p)>smoothing    % apply smoothing to traces
      p = smooth(p, smoothing, smoothMethod); 
    end
    else   % t1 empty or t2 empty?
      warning('trial %i had no timestamps');
      p=nan; t1=nan;
    end 
  else       % pos was empty. Indicates a dud trial.
    p=nan;   
  end
      if ~exist('scStartT','var') || isempty(scStartT) % no saccades on this trial?
        scStartT=nan; scVec=nan; scEndPt=nan; scBlink=true;
        scAmpl=nan; scBendT=nan; scBendR=nan; scVec=nan;
        scCurvA=nan; scCurvS=nan; scDepA=nan; scSpd=nan; scEndT = nan;
      end
        % now add all the saccades on this trial to the global list of saccades
      info.sRT    = nancat(1, info.sRT,    scStartT'-t1);
      info.sEndpt = nancat(1, info.sEndpt, scEndPt);
      info.sBlink = nancat(1, info.sBlink, scBlink);
      info.sAmpl  = nancat(1, info.sAmpl,  scAmpl);
      info.sBendT = nancat(1, info.sBendT, scBendT);
      info.sBendA = nancat(1, info.sBendA, scBendR); % was scBendA for angle. now use scBendR for distance
      info.sVec   = nancat(1, info.sVec,   scVec);
      info.sCurvA = nancat(1, info.sCurvA, scCurvA);
      info.sCurvS = nancat(1, info.sCurvS, scCurvS);
      info.sDepA  = nancat(1, info.sDepA,  scDepA);
      info.sSpd   = nancat(1, info.sSpd,   scSpd);
      info.sDur   = nancat(1, info.sDur,   scEndT' - scStartT'); % duration 
  
  %%%% Browsing saccades
  % plot traces trial-by-trial and then pause - very useful for debugging.
  if browse 
    if ~pupil && length(scStart)>0
      subplot(3,1,1);plot(abs(vs));
      ylim([0,100]);title(sprintf('trial %d',i));
      subplot(3,1,2);plot(p,'.'); xlim([-300 300]);ylim([-300,300]);
      hold on;
      plot(real([p(scStart) p(scEnd)])', imag([p(scStart) p(scEnd)])','b:');
      vecs=[p(scCurvT+scStart.'), p(scCurvT+scStart.')+10*scCurvS.'];
      plot(real(vecs)', imag(vecs)' , '-r' ); hold off;
      vecs=[p(scBendT+scStart.'), p(scBendT+scStart.')+40*scBendA.']; hold on;
      plot(real(vecs)', imag(vecs)' , '-g' ); hold off;
      subplot(3,1,3);plot([real(p), imag(p)]);legend({'X','Y'});
      ylim([-300,300]);
      hold on
      for(j=1:length(scStart))
        plot(scStart(j)*[1 1],get(gca,'ylim'),'b-');
        plot(scEnd(j)*[1 1],get(gca,'ylim'),'y-');
        plot( scStart(j)+[0 real(scVec(j))], [0 imag(scVec(j))], 'r');
        plot( [1 1]*(scStart(j)+scCurvT(j)),  [0 scCurvA(j)], 'g');
      end
      hold off;
    end
    if pupil
      plot([real(p) real(original_p)]);
    end
    title(sprintf('t %g',i)); 
    if     browse == 3,  keyboard % use this to stop and allow inspection
    elseif browse == 2,  pause;   % use this to allow keypress-to-advance
    end
  end
  
  %%%% Return requested samples
  if RETURN_TRACES % Do we want to return the individual traces for trials?
    if isempty(p) p=nan;end;
    % pack into results array
    if(stretch==0)               % no stretch: add nans at the end
      r=nanassign(r,[i,nan],p ); % R ( i, : ) = P.'    store trace in R
    else                         % stretch: interpolate data
      if(length(p)==1)           % so that each trial appears tp
        r(i,1:stretch)=p;                % have the same number of samples.
      else                       % if only one sample, duplicate it.
        %r(i,:)=interp1( s(i).pos(filter,1)-t1, p, linspace(1,length(p),stretch)' );
        r(i,:)=interp1( 1:length(p), p, linspace(1,length(p),stretch)' ); % edited 2014 july
      end                        % otherwise use interpolation to stretch.
    end
    
    %%%% Baseline subtraction
    if isnumeric(baselinesub) && baselinesub>0  % find baseline
      if length(p)>=baselinesub               % provided we have enough samples
        baseline=nanmean(p(1:baselinesub)); % calculate mean over baseline period
      else
        baseline=nanmean(p);                % (otherwise use whatever we have)
      end
    elseif ischar(baselinesub)                  % The baseline can be a field
      t=s(i).pos(:,1)>=s(i).(baselinesub);    % in which case just use
      if ~isempty(t)                          % the sample at that moment in time
        baseline=p(t(i));                   % as the baseline.
      else % something has gone wrong - couldn't find specified time.
        baseline=p(1);                      % just use the first timepoint
        warning(['snip: cant find time ' baselinesub]);
      end                                     % of the trial, and warn.
    else baseline=0; % no baseline
    end
    r(i,:)=r(i,:)-baseline;                     % subtract baseline from all samples
    
  end
end
if size(r,1) < numel(s) % if last trials not processed due to errors, add nans
  r(size(r,1)+1:numel(s),:) = nan;
end
if verbose % inform user about cleaning?
    fprintf('cleanup: %d short data segments removed from %d trials\n', sum(removedsegments), sum(removedsegments>0));
end
% restore warnings
warning(prev_warn_state.state, prev_warn_state.identifier);
%%%% Means
if RETURN_TRACES & ~isempty(meancriteria)    % use criteria to take means?
    ucr = unique(meancriteria);              % unique levels of the criterion.
    for(j=1:length(ucr))                     % select trials with a certain criterion
        f=meancriteria==ucr(j);              % mean of eye data is taken for 
        if(excludenan) f=f & ~any(isnan(r),2)'; end; % exclude trials with nans?
        cm(j,:)=nanmean(r(f,:));             % each time point across selected trials
        cs(j,:)=nanstd(r(f,:))/sqrt(sum(f)); % standard error of mean
        leg{j}=num2str(ucr(j));              % create legend entries for criteria
        if verbose fprintf('criterion %s has %d valid trials\n', leg{j}, sum(f)); end;
        if smoothing
          cm(j,:) = smooth(cm(j,:), smoothing, 'moving');
          cs(j,:) = smooth(cs(j,:), smoothing, 'moving');
        end
    end
    r=cm;                             % return the mean traces to the user
      
    info.sem=cs;                      % standard error for criteria
    info.criteria=ucr;                % and the criterion levels
    %%%% Plotting
    if(doplot)                        % plot means of all meancriteria?
        held=ishold();
        plot(cm.','-');
        hold on
        plot((cm-cs).',':');          % error bars as dots
        plot((cm+cs).',':');
        legend(leg);
        if(~held) hold off;end;
    end
else                                  % not taking means?
    if(doplot)
        isheld=ishold();              % basic plot of all saccades
        if(excludenan)                % if excluding nans from the plot, 
            f=~hasnan;                % select which trials to exclude
            if verbose fprintf('plot: %d trials removed containing nan\n', sum(hasnan)); end;
        else  f=ones(size(r,1),1)==1; % if not, include all trials
        end
        if pupil     plot(real(r(f,:))'); 
        else         plot(r(f,:).');
        end
    end
end % end of return traces.




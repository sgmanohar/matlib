function [squeezeData keys timeOfLastAcquisition] = waitForForceData(ex,...
    timeOfLastAcquisition, ...
    maxTimeToWait, stopRecordingAtThreshold , ISI, ...
    continuousFeedbackFunction )
% Function to stream data from squeezy device
% require MP Acquisition daemon to be running
% 'timeOfLastAcquisition':
%  Needs to know the time of the last acquisition in order to calculate 
%  how many samples to dequeue before starting this acquisition.
%  The first time you call this in a trial, you should send the time that
%  the acquisition started - i.e. tr.startSqueezyAcquisition.
%  Subsequent calls in the same trial must pass the value returned from
%  this function 'timeOfLastAcquisition'.
%
% The function blocks execution until either 
% 1) the force exceeds 'stopRecordingAtThreshold'
% 2) the time elapsed exceeds 'maxTimeToWait'.
%
% The function returns a matrix 'sqeezeData(TIME,CHANNEL)' 
% which is the force data for the left and right channels
% and also if any keys were pressed, they are in the vector 'keys',
% otherwise this vector is empty.
%
% Finally, we also return "timeOfLastAcquisition" which is the time that we
% last checked the MP150 stream. This is needed for further calls, to know
% how much data to de-queue
%
% continuousFeedbackFunction( currentForce )
%  should be a function that can be called repeatedly each time new samples
%  are acquired. You can do the drawing in here.


EXIT = false;
startRecordingTime = GetSecs;

lastread=GetSecs;

if ex.useSqueezy

%%%% De-queue data from MP150
% first dequeing is for all samples since the start of the acquisition.



% Calculate how many samples need to be read from the MP150; sadly we
% have to use  (last read-time) * (sample-rate) due to awful programming
% on behalf of BioPAC.
queueLen = @(lastread)2*floor((GetSecs-lastread)*ex.MP_SAMPLE_RATE);
% set up recording buffer of 3 seconds, and buffer
record=nan*zeros(ex.MP_SAMPLE_RATE*ISI,2);
lq = queueLen(timeOfLastAcquisition);
buffer=nan*zeros(lq,1);
[z buffer nread] = calllib(ex.mplib, 'receiveMPData', buffer,lq,0);
%lastread=timeOfLastAcquisition+nread/2/ex.MP_SAMPLE_RATE; % update last read time.
record(1:(nread/2),:)=[buffer(1:2:nread) buffer(2:2:nread)]; % store both channels
nTotalRead=nread/2; 
fprintf('err=%g',GetSecs-lastread);

% not sure how to code that, or should I use button press and if they responded 'yes' then 'get ready to squeeze the left grip with your left hand'
while(GetSecs < startRecordingTime + maxTimeToWait && ~EXIT)              % how do I do this for right and left hand
    %%%% CHECK KEYPRESSES
    [keyisdown,secs,keycode] = KbCheck; % check for real key
    keys=find(keycode);
    if(any(keys==27)) EXIT=1; end         % check for ESCAPE
    
    %%%% POLL squeezes from MP150
    lq = queueLen(lastread); % required length of queue

    if(lq>0)                 % any data in queue?
        if(lq>length(buffer)) buffer=nan*zeros(lq,1); end % extend buffer if needed
        [z buffer nread]=calllib(ex.mplib, 'receiveMPData', buffer,lq,0); % READ MP150
        %lastread=lastread+nread/2/ex.MP_SAMPLE_RATE; % keep log of elapsed reads, to allow calculation of queue length next time.
        
        lastread=GetSecs;
        %times(nTotalRead)=GetSecs; % keep track of time
        record((nTotalRead+1):(nTotalRead+nread/2),:)=[buffer(1:2:nread) buffer(2:2:nread)]; % store both channels
        nTotalRead=nTotalRead+nread/2; % keep track of how many bytes read
        if nTotalRead>20
            mu  = mean( record((nTotalRead-20):nTotalRead, :) );       % 20ms mean force
            mud = mean( diff(record((nTotalRead-20):nTotalRead, :)) ); % 20ms mean differential
            if any(mu > stopRecordingAtThreshold)
                EXIT=1;
            end
        end
    end
    if exist('continuousFeedbackFunction','var')&& exist('mu','var')
        continuousFeedbackFunction( mu ); % if the user supplied a feedback 
        % function, call it now with the mean force over last 20ms (this is
        % a vector with two elements, one for each channel)
    end
end

elseif ex.useLabjack
  % reading from the labjack is considerably simpler than reading from the
  % MP150 - since we can get away without queueing I think?
  record = [];
  
  currtm = GetSecs;
  nTotalRead = 0;
  while EXIT == false
    % Execute the requests.
    ex.ljudObj.GoOne(ex.ljhandle);
    
    % Get all the results.  The input measurement results are stored. All
    % other results are for configuration or output requests so we are just
    % checking whether there was an error.
    xtm=GetSecs-currtm;
    [ljerror, ioType, chan1, v0, dummyInt, dummyDbl] = ex.ljudObj.GetFirstResult(ex.ljhandle, 0, 0, 0, 0, 0);
    [ljerror, ioType, chan2, v1, dummyInt, dummyDbl] = ex.ljudObj.GetNextResult( ex.ljhandle, 0, 0, 0, 0, 0);

    nTotalRead = nTotalRead + 1; % keep track of how many samples read
    record(nTotalRead,:) = [ v0, v1 xtm ]; % changed from xtm for consistency
    
    if xtm > maxTimeToWait
      EXIT = true; 
    end
    pause(0.002) % wait 2 milliseonds  = 500 Hz

    if nTotalRead>20
      mu  = mean( record((nTotalRead-20):nTotalRead, :) );       % 20ms mean force
      mud = mean( diff(record((nTotalRead-20):nTotalRead, :)) ); % 20ms mean differential
      if any(mu > stopRecordingAtThreshold)
        EXIT = true;
      end
    end
    
    [keyisdown,secs,keycode] = KbCheck; % check for real key
    keys=find(keycode);
    if(any(keys==27)) EXIT=1; end         % check for ESCAPE
    if keys(KbName('LeftArrow'))
      record(end+1,:) = [ inf 0 0 ]; % like a squeeze on channel 1
      EXIT = 1;
    end
    if keys(KbName('RightArrow'))
      records(end+1,:) = [ 0 inf 0 ]; % on channel 2
      EXIT = 1; 
    end
    
    
    if exist('continuousFeedbackFunction','var')&& exist('mu','var')
      continuousFeedbackFunction( mu ); % if the user supplied a feedback
      % function, call it now with the mean force over last 20ms (this is
      % a vector with two elements, one for each channel)
    end

  end  %  while requestedExit == false
  
  
else % no squeezy or labjack? 
  % keyboard dummy version (for testing)
  keys = 1;  % wait for all keyboard keys to be releaed
  while any(keys),
    [~,~,keys] = KbCheck();
  end
  record=[]; 
  t = GetSecs;
  % wait for a keypress
  while ~any(keys) && (GetSecs < t+maxTimeToWait)
    [~,~,keys] = KbCheck();
    if keys(KbName('LeftArrow'))
      record = [record; [inf 0 0] ];
    elseif keys(KbName('RightArrow'))
      record = [record; [0 inf 0] ];
    else
      record = [record; [0 0 0] ];
    end
    pause(0.002);
    if exist('continuousFeedbackFunction','var')
      continuousFeedbackFunction( [1 1] ); % if the user supplied a feedback
    end
  end
end



squeezeData = record; % send back the actual data recorded
timeOfLastAcquisition = lastread; % send back the time we last read the device





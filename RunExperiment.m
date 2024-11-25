function result = RunExperiment ( doTrial, ex, params, blockStart )
% result = RunExperiment ( @doTrial, ex, params [, @blockStart] )
% The body of the experiment.
% Sanjay Manohar 2008
%
% What It Does:
% ------------
% 1) Create screen (type 'help prepareScreen' for details)
%
% 2) Initialise eyelink (if ex.useEyelink==1), recording default eye
%    position information to a date-numbered file. (if useEyelink==2 then
%    use a dummy eyelink setup)
%
% 3) Combine information from 'ex' and 'params' structures. Params is 
%    either interpreted as values that override those of ex, or alternatively
%    if params is the result from a previous experiment, then simply continue
%    the previous experiment.
%    If not continuing a previous experiment, create trial structure 
%    (type help createTrials) using parameters in ex.
% 
% 4) Iteratively call the given trial function: doTrial(scr, el, ex, trial) 
%    Parameters sent: screen structure, eyelink structure, 
%      experiment parameters structure, and trial parameters structure
%
% 5) Handles eyelink trial, calibration, drift correction, and abort logic.
%
% 6) Handles keypresses 
%      'R'    repeat trial now
%      'D'    drift correct
%      'C'    tracker setup screen
%      'F1'   allow keyboard input to modify experiment
%      escape exit block and exit experiment. 
%    The key result should be passed
%    back from doTrial in the return structure as trial.key, and trial.R
%    should return a status (R_ERROR terminates block).
% 
% 7) Handles errors: saves all data in 'errordump.mat', cleans screen and
%    shuts down eyelink; also attempts to upload EDF file to current
%    directory. Also, data is saved at every block-end, in LastExperiment.dat
%
% NOTE:
%  1) ex.R_ERROR (any number) and ex.useEyelink (true/false) must be defined.
%  2) doTrial must be a function handle (passed using @) with the
%     parameters as described in 4 above
%
% Input arguments:        
%                   
% doTrial:          user-supplied function of the form 
%                   trialResult = doTrial(scr, el, ex, trial)
%                   Called with parameters of each trial.
% 
%
% ex:               Experimental parameters, containing fields:
%  blocks           = number of blocks, or alternatively an array of 
%                     integers specifying the block types for each block to
%                     run. Default = 1
%
%  blockLen         = number of trials per block. 
%
%  blockVariables.varName1  = [value1, value2, value3]
%  trialVariables.varName2  = [value4, value5]
%                     Specify trial types. 
%                     These create variables that are shuffled across
%                     trials and blocks. Each variable 'varName_i' will take 
%                     one of the specified value_i from the corresponding 
%                     array. 
%                     You can use a cell array instead.
%                     
%  useEyelink       = whether to use eyelink
%  useScreen        = whether to use the PsychToolbox Screen command
%  useSqueezy       = initialise squeezy device. default initialises to 
%                     500Hz, 2 channels
% 
%  exptStart        = @function exptStart(scr, el, ex) is called at
%                     the start of the experiment. Can be used for
%                     one-off instructions.
%
%  blockStart:      = user supplied function of the form 
%                     @function blockStart (scr, el, ex, trial), where the trial
%                     supplied is the first trial of the block. use this to
%                     write a display needed at the start of a block.
%
%  randomSeed       = integer > 0: set the trial randomisation seed
%                     if omitted or 'nan', then randomise from the timer.
%
%  rethrowErrors    = if true, then the errors will be shown as proper
%                     errors, after the cleanup
% 
%  practiceTrials   = number of practice trials. If this is a multiple of
%                     the number of trial types, there will be one of each
%                     type.
%                     If this is a vector, it will use those trials -
%                     please note that it will wrap the values across
%                     blocks e.g. if you have blocks of 10 trials, trial 12
%                     will be block 2 trial 2. The block order will be the
%                     default order (1:nBlocks). This can be used to set
%                     specific trials, if shufflePracticeTrials is off (see
%                     below).
%
%  shuffleTrials    = true if you want the trials in a random order in each 
%                     block. Otherwise they will run in the order that the
%                     fields of 'trialVariables' were created. 
%
%  shufflePracticeTrials = (same but for practice trials)
% 
%  repetitionsPerBlock: You can specify this instead of blockLen. If this
%                     is 1, then the block length will be the same as the
%                     number of trial types.
%
%  MP_SAMPLE_RATE   = sample rate for the squeezy - Hz. Default 500 Hz
%  
%  
%
% params:           Either - additional experimental parameters which 
%                   override those in 'ex', or 'result' returned from
%                   a previous RunExperiment.
%                   set a field 'overrideParameters' if you want to
%                   automatically override any parameters set in the
%                   experiment program, otherwise you will be asked.
%
% blockStart:       can also be specified as a separate parameter, as an
%                   @function.
%
% in doTrial, you have access to 
% 
%  ex.R_NEEDS_REPEATING  }
%  ex.R_ERROR            }  = constants you can return as results in tr.R
%  ex.R_ESCAPE           }
%  ex.exitkey       the platform-specific Escape key number, for KbCheck
%
% result:
%  file             = EDF file name for eyelink data (hopefully transferred
%                     into local working directory at end of experiment)
%  trials(block,trial) = trial parameters structure - the specific 
%                     parameters created by createTrials, and sent 
%                     sequentially to each doTrial.
%  params           = The experimental parameters, combined
%  data(block,trial)= The results of each trial, as returned by doTrial.
%                     Once you have discarded unwanted trials/reordered the
%                     trials for analysis, you can use
%                     transpIndex(result.data) to access the values as
%                     matrices.
%  practiceResult(i)= results returned from doTrial for the practice
%                     trials.
%  last (1x2 double)= last block and trial that was successfully completed;
%                     this is the point to continue from if 'result' is
%                     used as the 'params' input to RunTrials.
%  date             = date/time string of start of experiment.
%
% Sanjay Manohar 2008

% setup parameters

% values of responseType, supplied by doTrial.
% send these values in tr.R to indicate the outcome of each trial
ex.R_ERROR        = -99;            % trial was an error; leave it and do nothing
ex.R_NEEDS_REPEATING = -97;         % error: rerun trial immediately after
ex.R_NEEDS_REPEATING_LATER = -95;   % error: rerun trial at end of block
ex.R_ESCAPE       = -98;            % escape was pressed - exit immediately
ex.R_INCOMPLETE   = -96;            % trial didn't complete as expected
ex.R_UNSPECIFIED  = -94;            % experiment didn't provide a return value
% Store the stack trace - so we know which .m experiment file was run 
% to execute the experiment
ex.experimentStack= dbstack;        
% Store the file record for top-level function (i.e. experiment) 
%  - includes the modification date and size! (2018) 
ex.experimentFile = dir(which(ex.experimentStack(end).file)); 

last = [1 1];                       % block and trial to start at
% fields from inputs that should not override current params
ignoreParamFields = {'experimentStack','experimentFile'};
if exist('params','var')           % if user specified a set of experimental 
    p=params;                       % parameters in the input parameters,
    if isfield(p,'params'), p=p.params; end
    fnames=fieldnames(p);           % override parameters from the original (ex)
    for x = 1:length(fnames);       % with those from the input structure (p).
        % ignore fields?
        if any(strcmpi(fnames{x}, ignoreParamFields))
          continue
        end
        override = 1;  % use values from params by default
        if isfield(ex,fnames{x})    % Is the field already in ex?
          if ~equals(p.(fnames{x}),ex.(fnames{x})) % if the values are unequal,
            if(isfield(p,'overrideParameters')) % if 'overrideParameters' is set, then use that value
                override=p.overrideParameters; 
                if(override)        % overwrite all parameters if 'overrideParameters' is specified
                    warning(['Using "' fnames{x} '" = "' p.(fnames{x}) '" from passed parameters.']);
                else
                    warning(['Using "' fnames{x} '" = "' ex.(fnames{x}) '" from experiment program.']);
                end
            else                    % otherwise ask specifically for each parameter
                try % try to print values nicely
                    if isa(p.(fnames{x}), 'function_handle') && isa(ex.(fnames{x}), 'function_handle') % if a function
                        override=(input(['Override ex "' fnames{x} '": "' char(ex.(fnames{x})) ...
                                '" with params: "' char(p.(fnames{x})) '" (1/0) ?'] ))
                    elseif isa(p.(fnames{x}), 'numeric') || isa(p.(fnames{x}), 'logical') % if numbers
                        override=(input(['Override ex "' fnames{x} '": "' num2str(ex.(fnames{x})) ...
                                '" with params "' num2str(p.(fnames{x})) '" (1/0) ?'] ))
                    elseif isa(p.(fnames{x}), 'char') % if characters
                        override=(input(['Override ex "' fnames{x} '": "' ex.(fnames{x}) ...
                                '" with params "' p.(fnames{x}) '" (1/0) ?'] ))
                    else % try this
                        override=(input(['Override ex "' fnames{x} '": "' ex.(fnames{x}) ...
                            '" with params "' p.(fnames{x}) '" (1/0) ?'] ))
                    end
                catch % if that doesn't work because of class, just display them
                    disp(fnames{x})
                    disp(ex.(fnames{x}))
                    disp(p.(fnames{x}))
                    override=(input(['Override ex "' fnames{x} '" with params? (see values above):  (1/0) ?'] ))
                end                 
            end
          end
        end
        if override
            ex.(fnames{x}) = p.(fnames{x});  % store field in ex
        end
    end
    if isfield(params,'last')   last   =params.last; end    % go straight to the last-executed trial 
    if isfield(params,'data')   results=params.data; end    % keep results of old trials
    if isfield(params,'trials') trials =params.trials; end  % keep old trial structure and randomisation 
    if isfield(params,'edfFiles')                          
        result.edfFiles=params.edfFiles;  % keep a list of EDF files used in previous runs
    else
        result.edfFiles={};
    end
    if isfield(params,'startTimes')       % keep a list of times that each run begins
        result.startTimes=params.startTimes;
    else
        result.startTimes={};
    end
    usingoldresults=1;
else usingoldresults=0; end

if ~isfield(ex,'useEyelink') ex.useEyelink=0; end  % default No Eyelink
if ~isfield(ex,'useScreen')  ex.useScreen=1;  end  % default Yes Screen
if ~isfield(ex,'useSqueezy') ex.useSqueezy=0; end  % default No Squeezy 
if ex.useEyelink, ex.useScreen=1;end               % can't have eyelink without screen
if ~ex.useScreen,  scr=0;end                       % not using screen?
if ~ex.useEyelink,  el=0;end                       % not using eyelink?
if ~exist('result','var'), result.startTimes={};end 
result.startTimes   = {result.startTimes{:}, datestr(now,31)}; % store time of experiment start
if ~exist('blockStart','var') && isfield(ex,'blockStart'), blockStart= ex.blockStart; end
if ~exist('exptStart','var') && isfield(ex,'exptStart'),   exptStart = ex.exptStart; end
if ~isfield(ex,'blocks'), warning('Assuming 1 block only'); ex.blocks=1;end
if ~isfield(ex,'deviceNumber'), ex.deviceNumber = -1; end % default is check all devices

% added 2016 to cater for different random number generator in new Matlab    
% (rand seed has been deprecated now)
if isfield(ex,'randomSeed') && numel(ex.randomSeed)==1 && ~isnan(ex.randomSeed)
  try   % if random seed present, set the fixed seed 
    rng(ex.randomSeed);
  catch mx   % older versions of matlab:
    rand('seed',ex.randomSeed);
  end
else    % no seed present: randomise the generator from the clock
  try
    rng('shuffle');
  catch mx  % older versions of matlab: 
    rand('seed',sum(100*clock));
  end
end

%%%%%%%%%%% CREATE TRIALS %%%%%%%%%%%%%%%%%
if exist('trials')~=1,  trials = createTrials(ex); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e=[]; % this carries any errors
fatal_error   = 0;
try     
    % Enable unified mode of KbName,
    % so KbName accepts identical key names on all operating systems:
    KbName('UnifyKeyNames');
    if ~isfield(ex,'exitkey')
      ex.exitkey = KbName('Escape');
    end
    if ex.useScreen
        scr=prepareScreen(ex);         % initialise screen (scr struct)
        ex.screenSize=scr.ssz;  
    end;
    
    if ex.useSqueezy                   % initialise squeezy device
      ex.mplib = 'mpdev'; mpdir='.';
      if ~isfield(ex, 'MP_SAMPLE_RATE') ex.MP_SAMPLE_RATE=500;end
      x=loadlibrary([mpdir '/mpdev.dll'],[mpdir '/mpdev.h']);
      [retval, sn] = calllib(ex.mplib,'connectMPDev',101,11,'auto');
      if ~strcmp(retval,'MPSUCCESS')
        fprintf('for eyelink use IP 100.1.1.2; for MP150 use IP 169.254.111.111\n');
        error('could not initialise MP150: %s', retval);
      end
      calllib(ex.mplib, 'setSampleRate', 1000/ex.MP_SAMPLE_RATE); % ms between samples
      calllib(ex.mplib, 'setAcqChannels', int32([1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]));  % which of 16 channels?
    end
    
    
    if ex.useLabjack
      ex.ljasm = NET.addAssembly('LJUDDotNet');
      ex.ljudObj = LabJack.LabJackUD.LJUD;
      % Read and display the UD version.
      disp(['UD Driver Version = ' num2str(ex.ljudObj.GetDriverVersion())])
      % Open the first found LabJack U3.
      [ljerror, ex.ljhandle] = ex.ljudObj.OpenLabJackS('LJ_dtU3', 'LJ_ctUSB', '0', true, 0);
      % assignments are in the factory default condition.
      ex.ljudObj.ePutS(ex.ljhandle, 'LJ_ioPIN_CONFIGURATION_RESET', 0, 0, 0);
      
      % Request a single-ended reading from AIN0m  and AIN1.
      ex.ljudObj.AddRequestS(ex.ljhandle, 'LJ_ioGET_AIN', 0, 0, 0, 0);
      ex.ljudObj.AddRequestS(ex.ljhandle, 'LJ_ioGET_AIN', 1, 0, 0, 0);
      ex.LJE_NO_MORE_DATA_AVAILABLE = ex.ljudObj.StringToConstant('LJE_NO_MORE_DATA_AVAILABLE');      
    end
    
    datestring=datestr(now,30);
    namebase=datestring([5:8, 10:13]); % base filename = date
    if(ex.useEyelink)                  % initialise eyelink (el struct)
        if Eyelink('IsConnected'), Eyelink('Shutdown');end;
        if(ex.useEyelink==1)
            success=Eyelink('Initialize');
            if success<0, fprintf('ensure eyelink software is on, cable is connected, and that\nthe IP address is set to 100.1.1.2\n'); return; end
        else
            Eyelink('InitializeDummy');
        end                            % set up calibration, record & receive gaze + pupil
        Eyelink('command', 'calibration_type = HV9');
        Eyelink('command', 'saccade_velocity_threshold = 35');
        Eyelink('command', 'saccade_acceleration_threshold = 9500');
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
        Eyelink('command', 'file_sample_data = GAZE,AREA,STATUS');
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,BUTTON');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,AREA');
        Eyelink('command', 'screen_pixel_coords = 0,0,%d,%d', scr.ssz(1), scr.ssz(2));
        Eyelink('Command', 'pupil_size_diameter = YES');
        el=EyelinkInitDefaults(scr.w);
        el.file=[namebase '.edf'];     
        Eyelink('openfile', el.file ); % create the EDF file
        result.file   = el.file;       
        el.disallowEarlySaccades=1;    % store the filename in output
        el.backgroundcolour=ex.bgColour;
        el.foregroundcolour=ex.fgColour;
        el.callback = []; % disable the flashy calibration screen
        if(~isfield(result,'edfFiles')) result.edfFiles={};end;
        result.edfFiles={result.edfFiles{:}, result.file};
        EyelinkDoTrackerSetup(el);     % run the calibration routine
        FlushEvents;                   
    end;
    result.trials = trials;     % save trial structure in output
    result.params = ex;         % save experimental parameters in output

    
    % if we are sending synchronisation triggers to EEG and eye tracker, send
    % end code
    if isfield(ex, 'sendTriggers') && ex.sendTriggers && isfield(ex.trigValues, 'elKeyword')
      trFake.block = 0; trFake.trialIndex = 0; % these needed for logevent
      TeenseyTrigger(ex.trigValues.startExpt, ex.trigValues.startExpt,...
             trFake, ex, el, 'startExpt');
    end

    % call experiment start code, if provided
    if exist('exptStart')
        exptStart(scr,el,ex);
    end

    % Practice trials
    % if there are practice trials, and we're not continuing from before:
    if isfield(ex,'practiceTrials') && prod(last)==1 
        % create a new set of random trials for the practice, in the same way as
        % would be done for the real experiment. 
        ex_prac = ex;
        if(isfield(ex,'blockorder'))
          ex_prac = rmfield(ex_prac,'blockorder'); 
        end
        ex_prac.blocks = ex.blocks; 
        %ex_prac.blockLen = ex.practiceTrials; 
        if ex.practiceTrials ~= ex.blockLen
            warning('practice doesnt contain all trial types');
        end
        if isfield(ex, 'shufflePracticeTrials')
            ex_prac.shuffleTrials = ex.shufflePracticeTrials; 
        end
        ex_prac.blockorder = 1:ex.blocks;
        % run create trials using the 'blockLen' of the practice trials. 
        % if this is a multiple of the number of trial types, then there
        % will be one of each trial type.
        prac=createTrials(ex_prac); 
        
        if ex.useScreen % show "Practice" screen, and wait for keypress
            drawTextCentred(scr, 'Practice trials', ex.fgColour);
            Screen('Flip', scr.w);
        else
          fprintf('Practice trials\n'); 
        end
        KbWait(ex.deviceNumber); 
        disp(prac)
		
        if numel(ex.practiceTrials) > 1%if ex.practiceTrials is a vector of trials to run
            if all(size(ex.practiceTrials)> 1)%if it is a matrix
                warning('ex.practiceTrials is a matrix. Practice trial order will run firstly by row, then by column');
                NP = reshape(ex.practiceTrials,1,numel(ex.practiceTrials));%flatten, by row then column
            else
                NP = ex.practiceTrials ; % number of practice trials
            end
        else
            NP = 1:ex.practiceTrials; %from 1:number of practice trials
        end
        
        for i=1:numel(NP) % for each practice trial,
          % run practice trials in order from NP
          tr = prac(1+floor((NP(i)-1)/ex.blockLen),1+mod(NP(i)-1,ex.blockLen));
          tr.isPractice = true; % mark the trial as a practice.
          % Now run each practice trial, setting the block index as '0'.
          tr = runSingleTrialAndProcess(scr, el, ex, tr,doTrial,0,i); 
          if isfield(result,'practiceResult') % if there are already practice results, 
            % make sure the new practice trial structure is compatible
            [result.practiceResult,tr]=ensureStructsAssignable(result.practiceResult,tr);
          end
          result.practiceResult(i)=tr; % store result in "practiceResult" in order they were run
          if tr.R == ex.R_ERROR || tr.R == ex.R_ESCAPE % if there was an error, 
            fatal_error=1; break;                       % exit practice trials 
          end
        end
    end
            
            

    %%%%%%%%%%%%%%%%% Main blocks and trials %%%%%%%%%%%%%%%%%
    if ex.useScreen
        Screen(scr.w, 'FillRect',ex.bgColour, scr.sszrect);
        drawTextCentred(scr, 'Start of Experiment', ex.fgColour);
        Screen('Flip', scr.w);
    else
        disp('Start of experiment - press a key');
    end
    KbWait(ex.deviceNumber);               % press a key to start the experiment
    
    
    if length(ex.blocks)==1, bnum=last(1):ex.blocks;
    else bnum = ex.blocks(last(1):end);
    end                   % continue from last block
    for b=last(1):ex.blocks
        if exist('blockStart','var')    % call the blockStart method if supplied
            kcode = 1; while any(kcode) [z z kcode]=KbCheck(ex.deviceNumber); end;
            FlushEvents '';
            tr=trials(b,1); tr.block = b; % allow the block start code to know which block we are in
            blockStart(scr,el,ex,tr);
        end
        repeatLater = []; % trials to repeat at end of block
        for t=1:ex.blockLen
            if b==last(1) && t<last(2), continue; end; % skip through trials already done
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run single trial
            tr=runSingleTrialAndProcess(scr,el,ex,trials(b,t),doTrial,b,t);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist('results','var')  % to append previous struct, ensure all the trials have the same data fields
                [results, tr]=ensureStructsAssignable(results,tr);
            else results=[]; % (first trial)
            end
            results = [results tr];  % then append
            if tr.R==ex.R_ERROR || tr.R==ex.R_ESCAPE
              fatal_error=1; break; 
            end
            if(isfield(ex,'R_NEEDS_REPEATING_LATER') && tr.R==ex.R_NEEDS_REPEATING_LATER)
              repeatLater=[repeatLater, t]; 
              fprintf('Blk %d Trial %d will be repeated at end\n',b,t);
            end;
            result.data   = results;      % store the results on the output
            result.last   = [b,t];        % keep track of what the last complete trial is
            % write the data to disc in a temporary file after every trial!
            % this allows data to be recovered in from this file after a fatal crash
            save 'LastExperiment' result  
        end  %end of block
        % repeat-later trials... keep going until they don't need repeating.
        % add them to the end of the data, but put in the appropriate trial
        % index for where it would have been in the sequence.
        t=1;
        while (~fatal_error) && (t<=length(repeatLater)) 
            tr=runSingleTrialAndProcess(scr,el,ex,trials(b,repeatLater(t)),doTrial,b,repeatLater(t));
            tr.trialIndex=repeatLater(t);   % make the trial index the same as it should have been
            tr.isRepeated=1;                % new flag to signify repeated trials
            [results tr]=ensureStructsAssignable(results,tr);
            results = [results tr];
            if(tr.R==ex.R_ERROR | tr.R==ex.R_ESCAPE) fatal_error=1; break; end;
            % allow the trial to repeated more than once
            if(tr.R==ex.R_NEEDS_REPEATING_LATER) repeatLater=[repeatLater, repeatLater(t)]; end; 
            result.data   = results;
            save 'LastExperiment' result
            t=t+1;
        end
        % at the end each block (or when quit), save all trials with the full filename
        save([namebase '.mat'], 'result'); 
        if ex.useScreen && ~exist('blockStart','var')
            drawTextCentred(scr, 'End of block', ex.fgColour);
            Screen('Flip', scr.w);
            KbWait(ex.deviceNumber);                             % wait for keypress after each block
        end
        [z z kcode]=KbCheck(ex.deviceNumber);
        if kcode(ex.exitkey) || fatal_error,  break;  end
    end
    %%%%%%%%%%%%%  end of experiment %%%%%%%%%%%%%%%%%
catch exc                  % in case of an error  
  % handle labjack errors:
  if(isa(exc, 'NET.NetException'))
    eNet = exc.ExceptionObject;
    if(isa(eNet, 'LabJack.LabJackUD.LabJackUDException'))
      disp(['UD Error: ' char(eNet.ToString())])
    else
      disp(['.NET Error: ' char(eNet.ToString())])
    end
  end
  disp(getReport(exc))
  
  e=lasterror;            % display error message
  fprintf('Error : %s\n',...
    e.message);
  for(i=1:length(e.stack))
    disp(e.stack(i));
    save 'errordump';
    if exist('results','var')
      result.data=results; % and still give back the data so far
    end
  end
end

% if we are sending synchronisation triggers to EEG and eye tracker, send
% end code
if isfield(ex, 'sendTriggers') && ex.sendTriggers && isfield(ex.trigValues, 'elKeyword') 
  trFake.block = 0; trFake.trialIndex = 0; % these are needed for logEvent
  TeenseyTrigger(ex.trigValues.endExpt, ex.trigValues.endExpt,...
      trFake, ex, el, 'endExpt');
end

if(ex.useScreen)        
  Screen closeall;      % restore screen
end                     
if ex.useEyelink        % close eye data file then transfer
    if(Eyelink('isConnected'))
        Eyelink('closefile');  
        fprintf('Downloading edf...');
        if Eyelink('receivefile',el.file,el.file) < 0
            fprintf('Error in receiveing file!\n'); 
        end
        fprintf('Done\n');
        Eyelink('shutdown');
    end
end
if ex.useSqueezy        % close squeezy device
  calllib(ex.mplib, 'disconnectMPDev');
  unloadlibrary('mpdev');
  fprintf('disconnecting MP150\n');
end
if ex.useLabjack        % close labjack device
  ex.ljudObj.Close();
end


FlushEvents '';
ShowCursor;             % show the mouse cursor again

if fatal_error && exist('tr','var') && isfield(tr,'R') % if ending abruptly, explain why
  if tr.R == ex.R_ESCAPE            % was escape pressed?
    fprintf('Exiting -- Escape key pressed\n');
    fprintf('If this was unintentional, you might be able to resume the experiment by re-running it\n');    
    fprintf('and using the result (e.g. ans or result from LastExperiment.mat) as the input parameter.\n');    
  end
end

% rethrow errors? if enabled, this will require any user code after the
% experiment to catch the error if finalisation is required. this allows
% debugging directly into the location of the problem.
if isfield(ex,'rethrowErrors') && ex.rethrowErrors
  if ~isempty(e), rethrow(e); end
end


%%%%%%%%% EXIT HERE %%%%%%%%%



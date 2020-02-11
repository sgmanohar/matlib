function tr = TeenseyTrigger(trigValue, waitForReturn, tr, ex, el, eventName)
% TeenseyTrigger(trigValue, waitForReturn, tr, ex, el, eventName)
% Communicates with the serial port (send command, flush, read output)
% Sends a trigger to the serial port, a KEYWORD trigger to the eye link (if
% connected), and logs the eventName in tr.
% Inputs:
%           trigValue is the duration of the pulse to send to the Teensey 
%               device e.g. 100 will become 2po1w100m0w to send a 100ms 
%               duration pulse to pin 2 (pin = ex.pinNumber)
%           waitForReturn [optional] how many ms to wait (duration of
%               trigger) before returning from function - to prevent
%               overlapping triggers. [] = no wait
%           tr, ex, el, eventName are all from doTrial and if passed will
%               be used by LogEvent() to send a time code into tr. ex must
%               have ex.pinNumber, ex.trigValues.elKeyword


trigCode =  sprintf('%dpo1w%dm0w', ex.pinNumber, trigValue); % write the code

fprintf(ex.portHandle, trigCode); % send code

% send trigger to eye link
if(ex.useEyelink)
    Eyelink('message', '%s', [ex.trigValues.elKeyword ' ' num2str(trigValue)]);
    Eyelink('command', ['record_status_message ' ex.trigValues.elKeyword ' ' num2str(trigValue)]);
end

tr = LogEvent(ex, el, tr, eventName); % 
% tr.(eventName) = GetSecs(); % log time in tr

if isa(ex.portHandle, 'serial') % won't do this if just printing to console
    flushinput(ex.portHandle); % empty input buffer
    flushoutput(ex.portHandle); % empty output buffer
end

if exist('waitForReturn','var') && ~isempty(waitForReturn)
    WaitSecs(waitForReturn/1000); % wait until code has finished
end

end
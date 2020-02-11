function portHandle = TeenseyInit(portName)
% portHandle = TeenseyInit(portName)
% Opens a serial port session
% Input: portName is the name/address of the port e.g. 'COM3' for windows 
% or % '/dev/tty.usbmodem40886101' for Mac
% Output: portHandle will be the handle for that port
% 
% You can close the connection by fclose(portHandle);

% close any open communication interface objects
try
    fclose(instrfind) ;
end

% get the handle to the specified port
portHandle = serial(portName); 

% open that port
try
    fopen(portHandle);
catch ME % in case it didn't work the first time
    if (strcmp(ME.identifier,'MATLAB:serial:fopen:opfailed'))
        try 
            fclose(instrfind);
            fopen(portHandle);
        catch ME
            rethrow(ME);
        end
    end
end

end
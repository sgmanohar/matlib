function [lostFix, tr] = MonitorFixation(el, location, accuracy, duration, ex, tr)
% Eyelink monitor fixation - return flag if goes outside accuracy radius
% [lostFix, tr] = MonitorFixation(el, location, accuracy, duration, ex, tr)
% Passive monitoring of eye fixation, returns flag if fixation is lost
% sends a trigger also
%
% uses global DCO - drift correct offset
% 
% returns lostFix = 1 if fixation went outside of radius

global DCO;
if(~exist('DCO','var') || length(DCO)~=2) DCO=[0 0]; end

s=Eyelink('newestfloatsample');
if ~isstruct(s) % something went wrong accessing the eyelink - disconnected?
  EyelinkIsConnected = Eyelink('isconnected') % display this
  warning('waitfix:notconnected','WaitForFixation could not access the eyelink');
  result=0; return
end

p0=[s.time s.gx(el.eye)+DCO(1) s.gy(el.eye)+DCO(2)];
p = repmat(p0,5,1);

count = 0;
totaltime = 0;
lostFix = 0;

while totaltime < duration
  if(Eyelink('NewFloatSampleAvailable'))
    s=Eyelink('NewestFloatSample');
    if(isstruct(s))
      dt=s.time-p(end,1);
      if(dt>0) % if it's a later timepoint, add to p-list
        p = [p(2:end,:); s.time s.gx(el.eye)+DCO(1), s.gy(el.eye)+DCO(2)];
        count=count+1;
        if(count<=5) continue; end % take at least 5 samples
        totaltime=totaltime+p(end,1)-p(end-1,1);
      end
    end
  end
  if ~lostFix && norm(p(end,2:3)-location)>accuracy
      tr = TeenseyTrigger(ex.trigValues.lostFix, [], tr, ex, el, 'lostFix'); % send trigger about direction but don't wait for return    
      lostFix = 1;
  end
  
end
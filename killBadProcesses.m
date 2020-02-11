function killBadProcesses(ADMIN_PW)
% function killBadProcesses(ADMIN_PW)
% Kill a large number of Windows processes that might interfere with the
% smooth running of an experiment. By default, the network services,
% printing, encryption, updates, and firewall are turned off -- so make
% sure the network cable is disconnected.
% 
% To do the job well, we need the Administrator password for the system
% (i.e. for username "Administrator"). If this is not provided, we might not
% have privileges to kill tasks.
% 
% sgm 2015


%%%%%%%%%%%%%%%%%%
% List of bad processes
processes = {
  'Search*' % windows search indexing
  'unsecapp*'
  'IAStor*'
  'DCPSysMgr*'
  'FIH32*'
  'fsav*' % antivirus
  'fsdfwd*'
  'fsgk*'
  'FSHDLL*'
  'FSM*'
  'fsorsp*'
  'fssm*'
  'nv*'  % nvidia services - for changing graphics modes
  'sql*' % microsoft sql databases
  'smax4pnp.exe'
  'TdmNotify.exe'
  'TdmService.exe'
  'WLID*' % "windows live"
  'wmpnetwk.exe' % windows media player
  'SMSvcHost.exe' % TCPIP port sharing for .Net
  'PdfPro7Hook.exe'
  'PDVDDXSrv.exe' % power DVD
  'Bcm*' % broadcom / Dell secuirity Device and task 
  'wuauclt.exe' % windows update
  'jusched*' % java updates
  'SeaPort.exe' % Microsoft Bing bar
  % 'tv_*'; 'TeamViewer*'   % team viewer daemon
  'splwow64.exe' % print spooler
  % 'WmiPrvSE.exe' % windows management interface provider
  
  %%%% windows 7
  'swi_service.exe'
  'ALMon.exe'
  'swc_service.exe'
  'ALsvc.exe'
  'savservice.exe'
  'SAVAdminService.exe'
  'Management Agent NT.exe'
  'jhi_service.exe'
  'LMS.exe'

  'PrivacyIconClient.exe'
  'unsecapp.exe' % WMI
  
  'RouterNT.exe'
  'armsvc.exe' % adobe updater
  
  %%%%
  % Graphics  - AMD card
  % 'CalibrationTool.exe' % for touch?
  % 'CLIStart.exe'; 'CCC.exe'; 'MOM.exe' % catalyst
  % 'atieclxx.exe'; 'atiesrxx.exe' % AMD 
  'igfxpers' % persistence of video settings
  'igfxsrvc'
  'igfxtray'
  
  'DellOSDService*' % on screen display
  'hkcmd.exe' % intel hotkey commands
  'IAStorDataMgrSvc*' % configure hard disks
  'IAStorIcon*'
  'o2flash.exe' % flash card storage
  'InputPersonalization.exe' % handwriting recognition
  'iprnt*' % print services
  'spoolsv.exe'
  'IPROSetMonitor*' % motherboard monitor
  'MediaButtons.exe'
  'SftService.exe' % dell backup
  'TabTip.exe' % tablet input panel
  'Toaster.exe' % recovery
  

  % 'taskeng.exe'; 'taskhost.exe'
  % 'wisptis.exe' % touch services
  % 'WavesSvc64.exe' % Maxx audio
  'TrustedInstaller.exe'
  
  %%% net
  'xtsvcmgr'
  'micasad.exe'
  'Zen*'
  'nw*'
  'nzrWinVNC*'
  'ZES*' % Zen endpoint security agent
  };

services = {
  'wuauserv' % windows update service
  'WSearch' % windows search indexing
  'W32Time' % keeps clock updated
  'UmRdpService' % remote desktop
  'TermService'  % "
  'SessionEnv'   % "
  'WinRM'   % remote management
  'Spooler' % print
  'SENS' % event notification
  'winmgmt' % stop firewall, SMS agent, and WMI console
  % 'secologon' % secondary logon
  % 'gpsvc' % group policy
  % 'lmhosts' % tcpip
  % 'DPS' % diagnostic policy service
  % 'wmiApSrv' % "WMI performance adapter"
  'Schedule' % task schedulaer
  'Power' % power up/down service
  'PolicyAgent'
  'Mcx2Svc'       % windows media center
  'WMPNetworkSvc' % "
  
  % network stuff...
  'nsi' % network store interface
  'NlaSvc' % network location
  'netprofm' 
  'Netman'
  'LanmanWorkstation'
  'LanmanServer'
  'iphlpsvc' 
  'Dnscache'
  'Dhcp'
  'Browser' 
  'BITS' % "background intelligent transfer"
  'WwanSvc' % wireless
  'Wlansvc' % "
  'wlidsvc' % windows live
  'WinHttpAutoProxySvc'
  'WebClient'
  'wcncsvs' % "windows connect now"
  % 'Netlogon'
  
  'WinDefend' % firewall
  'MpsSvc'    % "

  'WPCSvc' % parental controls
  'WerSvc'        % error reporting
  'wercplsupport' % "
  'eventlog' 
  'Wecsvc' % event collection
  'wbengine' % windows backup 
  'WatAdminSvc' % windows activation
  'TapiSrv' % telephony
  
  'CscService' % offline files
  % 'CryptSvc' % cryptographics
  'CertPropSvc' % certificate propagation
  'AeLookupSvc' % application experience
  'wscsvc'   % windows security service
  'VaultSvc' % credentials
  
  
  };

% if the parameter 'admin_pw' is provided, then use administrator account
% to kill tasks
USE_ADMIN = exist('ADMIN_PW','var'); 

if USE_ADMIN 
  command = 'runas /noprofile /user:Administrator "taskkill /f /im %s"';
else
  command = 'taskkill.exe /f /im %s';
end
done=false;


% KILL ALL PROCESSES
while ~done % loop until all processes have failed to be killed at least once
  done=true;
  for i=1:length(processes)
    if ~USE_ADMIN
      [status,result]=system( sprintf(command, processes{i}) );
      if status==0
        fprintf('%s killed\n',processes{i});
        done=false; % rerun the cycle
      end
      result
    else
      [p,exitcode,msg]=syscmd( 'taskkill.exe', sprintf('/f /im %s',  processes{i}), ADMIN_PW );
      disp(msg);
      if exitcode==0, 
        fprintf('killed %s\n',processes{i});
        done=false;
      end
      
    end    
  end
end

% KILL SERVICES
for i=1:length(services)
  syscmd('net', sprintf('stop %s', services{i}),ADMIN_PW);
end

% disable windows updates
syscmd('reg','add "HKEY_LOCAL_MACHINE\SOFTWARE\Microsoft\Windows\CurrentVersion\WindowsUpdate\Auto Update" /v AUOptions /t REG_DWORD /d 1 /f ',ADMIN_PW);
syscmd('reg','add HKLM\Software\Microsoft\Windows\CurrentVersion\Policies\Explorer /v HideSCAHealth /t REG_DWORD /d 0x1', ADMIN_PW);

function outputfn(s,e)
% collect process output in the global variable lastExitMsg
global lastExitMsg
if ~isempty(e.Data)
  lastExitMsg = e.Data;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, exitcode, stdout]=syscmd(proc, args, pwtext, stdin )
% call a system command with administrator privileges!
global lastExitMsg 
DEBUG=1;
  if DEBUG
    fprintf('%s %s',proc,args);
  end
  if exist('pwtext','var')
      p = System.Diagnostics.Process;
      p.StartInfo.FileName = proc;
      p.StartInfo.Arguments = args;
      p.StartInfo.UserName = 'Administrator';
      pw=System.Security.SecureString();
      for j=1:length(pwtext); pw.AppendChar(pwtext(j)); end
      p.StartInfo.Password =pw;
      p.StartInfo.LoadUserProfile = false;
      p.EnableRaisingEvents = true;
      p.StartInfo.CreateNoWindow = true;
      p.StartInfo.UseShellExecute = false;
      p.StartInfo.RedirectStandardOutput = true;
      p.addlistener('OutputDataReceived',@outputfn);
      p.StartInfo.RedirectStandardInput =true;
      p.Start();
      sw = p.StandardInput;
      p.BeginOutputReadLine();
      if exist('stdin','var') % output to stdin of process
        sw.WriteLine(stdin);
      end
      sw.Close();
      p.WaitForExit();
      exitcode=p.ExitCode;
      p.Close();
      stdout=lastExitMsg;
  else % fall back on non-administrator command
    [exitcode stdout]=system(sprintf('%s %s',proc,args));
    p=[];
  end
  if DEBUG,
    fprintf('%s\n',stdout);
  end
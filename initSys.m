function sys = initSys()
    sys = struct(...
        'vers',version, ...   %   Matlab version (e.g. 5.2, 6.5 or 7.0)
        'dbflag',0, ... %   flag used to debug
        'lineWidth',2,...
        'markSize',4,...
        'grey',[.6 .6 .6],...
        'pltCol','bgrky',...
        'pltTyp','ovdsp',...
        'debugMon',0);  % display output during monitoring?

     sys.vers=str2double(sys.vers(1:3));  % Version specific
     if (sys.vers>6.0)
         warning off MATLAB:divideByZero
     else
         warning off
     end    
    
%     sys.pltCol   = 'bgrky';                        % Plot item color (bgrky)
%     sys.pltTyp   = 'ovdsp';                        %   marker (circ,trian,diam,squ)
%     sys.lineWidth= 2;                              %   lines
%     sys.markSize = 4;                              %   markers
%     sys.grey  = [.6 .6 .6];                        %   color
%     sys.dbflag= 0;                                 % Used for debbuging
%     sys.vers=version; 

    
    
    % GUI Variables
% windowTypes windowTyp dataTypes dataTyp windowHndl cmdStr isRunning ...
% pltCol pltTyp lineWidth markSize grey

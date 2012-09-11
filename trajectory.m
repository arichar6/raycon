function varargout = trajectory(t,y,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  TRAJECTORY -- Ray trajectories
%             Part of the RAYcON package
%
%  ToDo    Optimize computing
%           - try to call dispersion only once per step in 'traj'
%
%  A. JAUN, Alfven Laboratory, KTH, 100 44 Stockholm, Sweden
%  A.N. KAUFMAN, Lawrence Berkeley Laboratory, Berkeley, CA 94720, USA
%  E.R. TRACY, College of William & Mary, Williamsburg, VA 23187-8795, USA
%
%  (C) Version 7.0,  14-Aug-2006. All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global plasma rays
%
% ----- Routine control
switch flag
case ''                                 % Return dy/dt = f(t,y)
  varargout{1} = traj(t,y);
case 'init'                             % Return IC [timespan,y,options]
  if (strcmp(plasma.PROBL,'tok'))
    [varargout{1:3}] = inittok(t,y);
  elseif strcmp(plasma.PROBL,'oth')
    disp('trajectory: undefined PROBL'); pause
  end
case 'events'                           % Return [value,isterminal,direction].
   [varargout{1:3}] = events(t,y);
otherwise
  error(['Unknown flag ''' flag '''.']);
end


function dydt = traj(t,y)
% ===== Equations of motion ==================================================
%global PROBL MODEL TYPE NRAY odeDim stp MON kray0 r0 c omega
global rays plasma sys cnst
%
if (strcmp(plasma.PROBL,'tok'))                    % ----- Evolution of ODE
 if (strcmp(rays.TYPE,'Trj')||strcmp(rays.TYPE,'Con'))
  dydt=disp_eig(y,'Trj');
 else
  %dydt=dispertok(0.,y,y,y,rays.poldot,'Amp');
 end
elseif strcmp(plasma.PROBL,'oth')                  % Steve define your own
 disp('trajectory: undefined PROBL'); pause
end

if (t>rays.ht(1)+1E-3*(rays.ht(1)-rays.ht(2)) && ...        % ----- Evolution of monitors
    (strcmp(rays.TYPE,'Con')||strcmp(rays.TYPE,'Amp')))
  for n=5:-1:1
    rays.ht(n+1)=rays.ht(n);      % Time              % Shift New -> Old
    rays.hy(:,n+1)=rays.hy(:,n);  % Position
    rays.hp(:,n+1)=rays.hp(:,n);  % Polarization
  end;
  % Put new values into history vectors
  rays.ht(1)=t; 
  rays.hy(:,1)=y;                  
  pol=disp_eig(y,'Pol');
  rays.hp(:,1)=reshape(pol,1,numel(pol));
  
  % if there's enough data, calculate the derivatives
  if (rays.ht(6)>0)                       
   ip1=1;% zp1=rays.hy(:,ip1); pp1=rays.hp(:,ip1);
   i0 =2; z0 =rays.hy(:,i0);%  p0 =rays.hp(:,i0);  
   dtp=rays.ht(ip1)-rays.ht(i0 );
   im1=3;% zm1=rays.hy(:,im1); pm1=rays.hp(:,im1); 
   dtm=rays.ht(i0 )-rays.ht(im1);
%    zdp=zp1-z0; zdm=z0-zm1;
%    pdp=pp1-p0; pdm=p0-pm1;
   if (dtp>0 && dtm>0)
       z=z0;
      rays.stp=rays.stp+1;                         % Accept new level
%     zdot=(zdp+zdm)'/(dtp+dtm);
%     zddot=(2*(zdp/dtp -zdm/dtm)/(dtp+dtm))';
%     rays.poldot=(pdp+pdm)'/(dtp+dtm);
    if (strcmp(plasma.PROBL,'tok'))                % Compute monitors
      out=disp_eig(z,'Mon');
    elseif strcmp(plasma.PROBL,'oth')              % Steve define your own
      disp('trajectory: undefined PROBL'); pause
    end
    rays.MON(rays.stp,:)=cat(2,t,out);                % Store monitors
   end;
  end;
end;


function [value,isterminal,direction] = events(t,y)
% ===== Monitoring conditions =============================================
global rays plasma sys cnst
%
%if (rays.stp<6)                            % ----- Evolution of monitors

% give the ray a few steps to get going.
% if this is a converted ray, it needs to get out of the conversion region
if (rays.stp<15)                            % ----- Evolution of monitors
  rays.adjcnv    = zeros(rays.NRAY,1);      % No adjustments to start
  rays.adjinit=1;    
  rays.adjctc    = zeros(rays.NRAY,1);
  value     =+1;                            % Arbitrary sign to begin
  direction = 0;
  isterminal = 1;
  %cnv=1; ctc=1;
  rays.moncnv=1; rays.monctc=1; 
else
  pts=6;                                    % ----- Nbr points to fit
  all=size(rays.MON,1); z0=[]; z1=z0; z2=z0;     % initialize
  mon=(rays.MON(all-pts+1:all,2:end));           % mon = monitors(time)
  tim=rays.MON(all-pts+1:all,1);                 % time
  sca=tim(1); tim=tim/sca; tt=t/sca;        % rescale
  A=[ones(size(tim)) tim tim.^2];           % quadratic fitting matrix
  nmon=size(mon,2)/rays.NRAY;
  for kray=1:rays.NRAY
    os=(kray-1)*nmon;
    lastwarn('');
    s=warning('off','MATLAB:rankDeficientMatrix');
    cf=A\mon(:,os+1:os+nmon);               % coefficients of parabola
    warning(s);
    [msgstr, msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:rankDeficientMatrix')
        % this is something of a hack... all the recently calculated ray is
        % lost.
        rays.monErrAbort = 1; 
        error('trajectory:events',...
            'Error calculating conversion and caustic monitors. Abort ray tracing.')
    end
    z0=[z0;cf(1,:)+cf(2,:)*tt+cf(3,:)*tt^2];%   columns = monitors
    z1=[z1;2*tt*cf(3,:)+cf(2,:)];           %   lines = rays
    z2=[z2;2*cf(3,:)];
  end
  jcaust=1; rays.monctc=z0(:,jcaust);            % Zero value ==> caustic
  jconvt=2; rays.moncnv=z1(:,jconvt);            % Zero derivative ==> convert
%
%
  if (rays.adjinit)                         % ----- Initialize normalization
    for kray=1:rays.NRAY
      %os=(kray-1)*nmon;                     % Conversion
      if (rays.adjcnv(kray)>0 && abs(rays.adjcnv(kray))<1E-6)
        rays.adjcnv(kray)=-rays.adjcnv(kray);   %   change old sign as cross 0
      else
        rays.adjcnv(kray)=1./rays.moncnv(kray); %   update with current value
      end
                                             % Caustic
      if (rays.adjctc(kray)>0 && abs(rays.adjctc(kray))<1E-6)
        rays.adjctc(kray)=-rays.adjctc(kray);   %   change old sign as cross 0
      else
        rays.adjctc(kray)=1./rays.monctc(kray); %   update with current value
      end
    end % kray
    rays.adjinit=0;
  end
  cnv = prod(rays.adjcnv.*rays.moncnv);
  ctc = prod(rays.adjctc.*rays.monctc);
 
  % this has a problem with 'Con' calculation...size(y)=[1 4]
  %ctc =  (~rays.inKspace(1)*2500 + rays.inKspace(1)*0.01)- norm(y(5:7));
  %ctc =  (~rays.inKspace(1)*1000 + rays.inKspace(1)*0.01)- norm(abs(y(5:7)));
  rays.monctc = ctc;
            % modification by Steve
  
% Another monitor to check that the value of the dispersion function is not
% too far from zero
  %dsp = abs(dispertok(0.,y,y,y,rays.poldot,'Dsp'))-1000;
  %    

  value(1)      = cnv;
  direction(1)  = 0;
  isterminal(1) = 1;
%  if (rays.odeDim~=4)
%    value(2)      = ctc;
%    direction(2)  = 0;
%    isterminal(2) = 1;
%  end
  if (rays.odeDim~=4)
      value(1)      = cnv*ctc;
      direction(1)  = 0;
      isterminal(1) = 1;
  end


  % ----- Display evolution
  if(sys.debugMon)
      t1=sprintf('%i time=%e ',rays.stp,t);
      t2=[' con=' sprintf('%0.1e ',(rays.adjcnv.*rays.moncnv)')];
      t3=[' cau=' sprintf('%0.1e ',(rays.adjctc.*rays.monctc)')];
      if (strcmp(rays.TYPE,'Con'))
        disp([t1 t2])
      elseif (strcmp(rays.TYPE,'Amp'))
        disp([t1 t2 t3])
      end
  end

end






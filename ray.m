function out = ray(action)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  RAYcON --  Ray tracing with mode conversion in a tokamak
%
%  A. JAUN, Numerical Analysis, KTH, 100 44 Stockholm, Sweden
%  A.N. KAUFMAN, Lawrence Berkeley Laboratory, Berkeley, CA 94720, USA
%  E.R. TRACY, College of William & Mary, Williamsburg, VA 23187-8795, USA
%
%  Documented under "http://www.nada.kth.se/~jaun"
%
%  To Do:
%   - Maslov: correct sign(det(focus)) -> sign(eig(focus))
%   - Why Eminus -> 0 at cyclotron resonance?
%
%  (C) Version 7.0,  14-Aug-2006. All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ----- Global variables documentation ---------------------------------------
global ...     % Plasma configuration
    plasma     % A struct which will contain the system config variables
%     EQ ...     %   equilibrium (either 'Solovev' or a filename)
%     MODEL ...  %   dispersion equation model (msw1x1, cld2x2, cld3x3, etc)
%     PROBL ...  %   problem (tok=tokamak, oth=Steve define your own)
%     acharge ...%   atomic charge of species (e.g. -1 for electron)
%     amass ...  %   atomic mass of species (e.g. 4 for helium)
%     b0 ...     %   magnetic field on axis
%     elong ...  %   plasma elongation
%     freq ...   %   antenna linear frequency
%     iaspr ...  %   plasma inverse aspect ratio
%     psin ...   %   magnetic flux at the plasma edge for normalization
%     q0 ...     %   plasma safety factor on axis
%     r0 ...     %   plasma major radius
%     n0 ...     %   densities on axis
%     na ...     %   densities profile factor in n=n0*(1-na*s^2)^nb
%     nb ...     %   densities profile factor in n=n0*(1-na*s^2)^nb
%     t0 ...     %   temperatures on axis
%     ta ...     %   temperatures profile factor in t=t0*(1-ta*s^2)^tb
%     tb ...     %   temperatures profile factor in t=t0*(1-ta*s^2)^tb
%     kant ...   %   antenna wave vector in (r,phi,z)
%     sant ...   %   antenna normalized radius
%     thant ...  %   antenna lower and upper poloidal angles
%     omega ...  %   antenna circular frequency
%     s ...      %   normalized radial mesh
%     theta      %   poloidal mesh
%     depo ...   %   radial deposition profiles per species
global ...     % Raytracing variables
    rays       % 
%     NRAY ...   %   number of rays to propagate
%     MON ...    %   unsorted monitors accumulated during the evolution
%     MONNorm ...%   normalizing factors
%     TYPE ...   %   calculation (Trj=trajectory,Con=conversion,Amp=amplitude)
%     aThres ... %   amplitude threshold below which rays are removed
%     convert ...%   list of indices for rays to convert
%     hp ...     %   5 steps history for time derivatives, value polarizations
%     ht ...     %   5 steps history for time derivatives, value of t
%     hy ...     %   5 steps history for time derivatives, value of y=(x,k,S,A,T)
%     inKspace...%   flag for evolution of each ray (0=config space, 1=k-space)
%     monctc ... %   last recorded caustics monitor
%     moncnv ... %   last recorded conversion monitor
%     plt ...    %   components of ODE to plot
%     stp ...    %   index to monitors along the ray
%     caustic... %   list of indices for rays to perform Maslov transform
%     odeDim ... %   dimension of ODE vector
%     odeOptions %   controlling ODE solver
%     NE         %   event ray number (which ray converts)
%     TR ...     %   solution times (rays' trajectories)
%     YR ...     %   solution components (rays' trajectories)
%     time ...   %   current time
%     timespan ... % duration of complete calculation
%     timeintv ... % duration until next stop
%     tspan ...  %   time interval tspan = [time time+timeintv]
%     y   ...    %   ODE components
%     kray0 ...  %   ray pencil initial wave vectors in (r,phi,z)
%     sray0 ...  %   ray pencil initial normalized radii
%     thray0 ... %   ray pencil initial poloidal angles
global ...   % Constants [MKSA]
    cnst     %
%     c ...      %   speed of light in vacuum [m/s]
%     e ...      %   electron charge [C]
%     eps0 ...   %   permittivity of vacuum [F/m]
%     mp         %   proton mass [kg]
global ...   % Physical variables
    sys      %
%    vers ...   %   Matlab version (e.g. 5.2, 6.5 or 7.0)
%     dbflag ... %   flag used to debug
%               % GUI Variables
%  windowTypes windowTyp dataTypes dataTyp windowHndl cmdStr isRunning ...
%  pltCol pltTyp lineWidth markSize grey
%

% Check if the global structs need to be initialized
if isempty(plasma)
    plasma = initPlasma;
end
if isempty(rays)
    rays = initRays;
end
if isempty(cnst)
    cnst = initCnst;
end
if isempty(sys)
    sys = initSys;
end

%
% ----- Initialize without calling argument (GUI control) --------------------
if nargin < 1,
  winSetup;                                  % Open window
  ray('preset');                             % Constants
  ray('evalText');                           % Read window input
  ray('auxval');                             % Auxiliary quantities
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Sequence programming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% ----- Terminate run --------------------------------------------------------
% elseif (isempty(rays.NRAY) || rays.NRAY<1) % I don't like this, it
% prevents preset and auxval from working before 'data' is run...
%                                            % Return without doing anything
%
% ----- Begin tracing --------------------------------------------------------
elseif strcmp(action,'begin');               % Start tracing
  ray('auxval'); ray('start');
  figure(1);clf; ray('propagate'); ray('boundary');
  figure(2);clf; ray('accumul');   ray('boundary');
  figure(3);clf; ray('history3');
  figure(4);clf; ray('history4');
%
elseif strcmp(action,'other');               % Steve define your own
  plasma.PROBL='oth';
%
%
% -----  specific sequences ------------------------------------------
elseif strcmp(action,'continue');
  %clear global MON TR YR; stp=0;             % Only to prevent false events
  rays.MON=[]; rays.TR=[]; rays.YR=[]; rays.stp=0;
  ray('remove_idle');                        % Rays without power
  if strcmp(plasma.PROBL,'tok')                     % Tokamak specific
    figure(1);        ray('propagate');
    figure(2);hold on;ray('accumul')
    figure(3);hold on;ray('history3');
    figure(4);hold on;ray('history4');
  elseif strcmp(plasma.PROBL,'oth')                 % Steve define your own
     disp('ray: undefined PROBL'); pause
  else
    disp('ray: undefined PROBL'); pause
  end
%
%
% ----- Caustics -------------------------------------------------------------
elseif (strcmp(action,'caustic')  && strcmp(plasma.TYPE,'Amp'))
  ray('caustic_which');                      % Determine rays to transform
  ray('caustic_')                            % Proceed with transformation
  ray('continue')                            % Resume evolution
%
elseif (strcmp(action,'caustic_') && strcmp(plasma.TYPE,'Amp'))
  if (numel(rays.caustic)>0)
    ray('caustic_list');                     % Proceed with transform
    %clear global MON TR YR; stp=0;           % Re-initilize history
    rays.MON=[]; rays.TR=[]; rays.YR=[]; rays.stp=0;
  end
%
%
% ----- Conversion -----------------------------------------------------------
elseif (strcmp(action,'convert')  && ...
        (strcmp(rays.TYPE,'Con')||strcmp(rays.TYPE,'Amp')))
  ray('convert_which');                      % Determine rays to convert
  ray('convert_')                            % Proceed with conversion
  ray('continue')                            % Resume evolution
%
elseif (strcmp(action,'convert_') && ...
       (strcmp(rays.TYPE,'Con')||strcmp(rays.TYPE,'Amp')))
  if (numel(rays.convert)>0)
    ray('convert_list');                     % Proceed with conversion
    %clear global MON TR YR; stp=0;           % Re-initilize history
    rays.MON=[]; rays.TR=[]; rays.YR=[]; rays.stp=0;
  end
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Elementary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% ----- Preset default values ------------------------------------------------
elseif strcmp(action,'preset');              % Constants
    
    cnst.c        = 2.9979E+8;                      %   speed of light in vacuum [m/s]
    cnst.e        = 1.6022E-19;                     %   electron charge [C]
    cnst.mp       = 1.6726E-27;                     %   proton mass [kg]
    cnst.eps0     = 8.8542E-12;                     %   permittivity of vacuum [F/m]

%
%
% ----- Auxiliary quantities -------------------------------------------------
elseif strcmp(action,'auxval')                   % Auxiliary quantities

  rays.aThres   = 0.05;                          % Amplitude threshold
  if strcmp(plasma.PROBL,'tok')                  %  Tokamak specific
      plasma.psin=0.5*plasma.b0/plasma.q0* ...
      plasma.elong*(plasma.r0*plasma.iaspr)^2;   %   magnetic flux at edge
      plasma.omega = 2*pi*plasma.freq;           %   circular frequency [rad/sec]
      rays.sray0=plasma.sant*ones(1,rays.NRAY);  %   rays IC in (s,theta)
          if rays.NRAY==1
            rays.thray0=plasma.thant(1);
          else
              rays.thray0=plasma.thant(1):...
                (plasma.thant(2)-plasma.thant(1))/(NRAY-1):plasma.thant(2);
          end;
      rays.kray0=ones(rays.NRAY,1)*plasma.kant;       % Antenna choice not necessarily 
                                                      %   compatible with dispersion
      plasma.depo=zeros(20,length(plasma.amass));     % Deposition profiles
%
 elseif strcmp(plasma.PROBL,'oth')                  % Steve define your own
     disp('ray: undefined PROBL'); pause
 end
%
%
% ----- Plot 1D equilibrium quantities in the mid-plane ----------------------
elseif strcmp(action,'equil1D')
 if strcmp(plasma.PROBL,'tok')                      % Tokamak specific
  rays.odeDim=4;
  plasma.s=-0.9998:.999/(plasma.NS-1):1.000; 
  plasma.s=abs(plasma.s);    %  flux surfaces
  plasma.theta=cat(1,pi*ones(plasma.NS-1,1),zeros(plasma.NS,1))'; %  polar angle
  [rho,r,z]=mapFlux(plasma.s,plasma.theta);                %  coordinates
  d=cat(2,r',z',ones(size(r))'*plasma.kant(1),ones(size(r))'*plasma.kant(3))';
  D=dispertok(0.,d,d,d,0,'Frq');
  h=plot(r,D'./1E6,r,plasma.freq/1E6*ones(size(r))); set(h,'LineWidth',sys.lineWidth)
  xlabel('major radius [m]'); 
  ylabel('hybrid- and cyclotron frequencies [MHz]');
 elseif strcmp(plasma.PROBL,'oth')                  % Steve define your own
     disp('ray: undefined PROBL'); pause
 end
%
%
% ----- Root from full dispersion relation closest to 'msw1x1' ---------------
elseif strcmp(action,'roots1D')
 if strcmp(plasma.PROBL,'tok')                      % Tokamak specific
  rays.odeDim=4;
  plasma.s=-.988:0.97/(plasma.NS-1):1.; 
  plasma.s=abs(plasma.s);          %  flux surfaces
  plasma.theta=cat(1,pi*ones(plasma.NS-1,1),zeros(plasma.NS,1))'; %  polar angle
  [rho,r,z]=mapFlux(plasma.s,plasma.theta);                %  coordinates
  %k0=rays.kray0(1,:); 
  %k0(1)=0;
  y=cat(2,r',z',zeros(size(r))'*plasma.kant(1),ones(size(r))'*plasma.kant(3));
  yv=reshape(y',1,size(y,1)*size(y,2))';
  T=dispertok(0,yv,zeros(size(yv)),zeros(size(yv)),0,'Msw');
  for kk=1:size(r,2)
    if T(kk)>0
      kms(kk)=sqrt(2)*omega/c.*sqrt(T(kk));
      kr=[0 0 1 0]';                         %  direction of adjustment
      krt(kk)=fzero('dispertok',kms(kk), ...
                   optimset('Display','off'), ...
                   y(kk,:),kr,zeros(size(kr)),0,'Dsp');
      if ( abs(krt(kk)) > 5*abs(kms(kk)) ) krt(kk)=0; end;
    else
      kms(kk)=sqrt(-T(kk));
      krt(kk)=0;
    end;
  end
  hold on;
  h=plot(r,kms,r,krt); set(h,'LineWidth',lineWidth); 
% legend('approx','exact');
  h=plot(r,krt,'ko'); set(h,'MarkerSize',markSize);
%                     set(h,'LineWidth',lineWidth); 
  xlabel('major radius [m]'); ylabel('wave vector  k_R [m^{-1}]');
  hold off;view(2)
 elseif strcmp(PROBL,'oth')                  % Steve define your own
     disp('ray: undefined PROBL'); pause
 end
%
%
% ----- Plot 2D mesh for position in (s,theta) -------------------------------
elseif strcmp(action,'mesh')
 if strcmp(plasma.PROBL,'tok')                      % Tokamak specific
  ns=11; nt=72;                              %  ds=0.1 and dt=5 deg
  s=0.0001:(.999/(ns-1)):1.;                 %  radial mesh
  theta=0:2*pi/nt:2*pi;                      %  poloidal mesh
  for kk=1:ns
    si=s(kk).*ones(size(theta));             %  flux surfaces
    [rho,r(kk,:),z(kk,:)]=mapFlux(si,theta); %    coordinates
  end;
  colormap(sys.grey)
  h=mesh(r,zeros(size(r)),z);
%  colormap('default')
  xlabel('major radius  R [m]');  zlabel('vertical direction  Z [m]');
  view(0,0); 
 elseif strcmp(plasma.PROBL,'oth')                  % Steve define your own
     disp('ray: undefined PROBL'); pause
 end
%
%
% ----- Plot domain boundary -------------------------------------------------
elseif strcmp(action,'boundary')
 if strcmp(plasma.PROBL,'tok')                      % Tokamak specific
  s=1.;                                      %  radial mesh
  theta=0:2*pi/plasma.NT:2*pi;                      %  poloidal mesh
  for kk=1:1
    si=s(kk).*ones(size(theta));             %  flux surfaces
    [rho,r(kk,:),z(kk,:)]=mapFlux(si,theta); %  coordinates
  end;
  hold on; h=plot3(r,zeros(size(r)),z,'k');
  set(h,'LineWidth',sys.lineWidth,'Color',sys.grey); hold off
  xlabel('major radius  R [m]');  zlabel('vertical direction  Z [m]');
  axis('equal');axis('tight');view(0,0);
 elseif strcmp(plasma.PROBL,'oth')                  % Steve define your own
     disp('ray: undefined PROBL'); pause
 end
%
%
% ----- Plot 2D equilibrium quantities and dispersion manifold ---------------
elseif (strcmp(action,'equil2D')||strcmp(action,'roots3D'))
 if strcmp(plasma.PROBL,'tok')                      % Tokamak specific
  rays.odeDim=4; sa=0.95;                         %  shave external Langmuir mode
  s=0.0001:0.99*sa/(plasma.NS-1):sa;                %  radial mesh
  theta=0:2*pi/(plasma.NT-1):2*pi;                  %  poloidal mesh
  for kk=1:plasma.NS
   if strcmp(plasma.EQ,'Solovev')                   %  Solovev equilibrium
    si=s(kk).*ones(size(theta));             %   flux surfaces
    [rho,ri(kk,:),zi(kk,:)]=mapFlux(si,theta);%   coordinates
   else                                      %  Numerical equilibrium
    sprintf('ERROR: gradShafr(%s) interface not implemented',plasma.EQ);
    %[rho,ri,zi]=gradShafr(EQ);
   end;
  end;
  if (strcmp(action,'equil2D'))              % --- 2D Projection to (R,Z)
    r =reshape(ri,numel(ri),1);
    z =reshape(zi,numel(zi),1);
    y =cat(2,r,z,1.5*ones(size(r))*plasma.kant(1),1.5.*ones(size(r))*plasma.kant(3));
    yv=reshape(y',1,size(y,1)*size(y,2));
    Dv=dispertok(0.,yv,yv,yv,0,'Dsp');
    Di=reshape(Dv',plasma.NS,plasma.NT);
    surf(ri,zi,Di); shading interp;  
    if (sys.vers~=6.5)                           % work around bug
       view(0,90); axis('equal');
    end
    xlabel('major radius [m]'); ylabel('vertical position [m]');
    title('Dispersion function or asinh(D)');
    colorbar
  elseif (strcmp(action,'roots3D'))          % --- 3D dispersion manifold
    k=linspace(-abs(plasma.kant(1)),abs(plasma.kant(1)),plasma.NK);
    R=zeros(size(ri,1),size(ri,2),size(k,2)); Z=R; K=R; Di=R;
    r =reshape(ri,numel(ri),1);
    z =reshape(zi,numel(zi),1);
    for i3=1:size(k,2)
       R(:,:,i3)=ri; Z(:,:,i3)=zi; K(:,:,i3)=ones(size(ri)).*k(i3);
       y =cat(2,r,z,ones(size(r))*k(i3),ones(size(r))*plasma.kant(3));
       yv=reshape(y',1,numel(y))';
       Dv=dispertok(0.,yv,yv,yv,0,'Dsp');
       D(:,:,i3)=reshape(Dv',plasma.NS,plasma.NT);
    end
    rtcol=patch(isosurface(R,K,Z,D,0.)); 
%   set(rtcol,'FaceColor',[.9 .9 .9],'EdgeColor',[.8 .8 .8]);
    set(rtcol,'FaceColor','blue','EdgeColor','red');
    xlabel('major radius [m]'); ylabel('wave vector k_R [m^{-1}]'); 
    zlabel('vertical position [m]');title('Dispersion manifold');view(20,25);
  end
 elseif strcmp(plasma.PROBL,'oth')                  % Steve define your own
     disp('ray: undefined PROBL'); pause
 end
%
%
% ----- Initial conditions for ray propagation -------------------------------
elseif strcmp(action,'start')
  %clear global MON MONNorm TR YR;            %  reset monitored quantities
  rays.MON=[]; rays.MONNorm=[]; rays.TR=[]; rays.YR=[];%  reset monitored quantities
  rays.time=0; rays.stp=0; rays.inKspace=zeros(1,rays.NRAY);
  rays.timeintv=rays.timespan;
  [tspan,y,odeOptions] = feval('trajectory',rays.timeintv,[],'init');
  rays.tspan = tspan;
  rays.y = y;
  rays.odeOptions = odeOptions;
  for kray=1:rays.NRAY
    os=(kray-1)*rays.odeDim;
    disp(['Ray#' sprintf('%i  init =',kray)...
                 sprintf(' %0.3g',y(os+1:os+rays.odeDim)) ])
  end
%
%
% ----- Remove idle rays -----------------------------------------------------
elseif strcmp(action,'remove_idle') && strcmp(rays.TYPE,'Amp')
  %rays.aThres=0.05;                               % Amplitude below which deleted
  logThres2=log(rays.aThres^2); idx=[]; 
  newn=rays.NRAY; newy=[]; newk=[];
  for kk=0:rays.NRAY-1
    off=kk*rays.odeDim;
    if (rays.y(off+8)<logThres2 && ~rays.inKspace(kk+1))
      idx=[idx sprintf(' %i,',kk+1)];
      newn=newn-1;
    else
      newy=[newy;rays.y(off+1:off+12)];
      newk=[newk rays.inKspace(kk+1)];
    end
  end
  if (~isempty(idx))
    tot=sprintf(' %i',newn);
    disp(['==> removing ray nbr' idx ' remains a total of ' tot]);
    rays.NRAY=newn; rays.y=newy; rays.inKspace=newk;
  end
  if rays.NRAY<1
    disp('==> All the energy has been deposited. Run terminated.');
  end
%
%
% ----- Propagate ------------------------------------------------------------
elseif strcmp(action,'propagate')
  figure(1); legend off;                     % --- Resume propagation
  %timeintv=rays.timespan-rays.time;
  timeintv = rays.timeintv;
  %rays.initialstep = 0.0001*timeintv;
  [tspan,y,odeOptions] = ...                 % Intermediate conditions
       feval('trajectory',timeintv,rays.y,'init');
  rotate3d on;  zoom off;  watchoff;
  % Trace the ray!
  try
%      [tr,yr]=feval('ode23s','trajectory',tspan,y,odeOptions);
      [tr,yr]=feval('ode45','trajectory',tspan,y,odeOptions);
      rays.time=tr(end); rays.y=yr(end,:)';      % Save new data
      rays.InK = rays.inKspace + zeros(size(tr));
      rays.TR=cat(1,rays.TR,tr);     rays.YR=cat(1,rays.YR,yr);
      if (strcmp(rays.TYPE,'Con')||strcmp(rays.TYPE,'Amp')) % --- Monitoring
        if strcmp(rays.TYPE,'Amp')
          ray('caustic_which');                  %   caustic candidates
        end
        ray('convert_which');                    %   conversion candidates
      end
      if ((tspan(2)-tr(end))/tspan(2)<0.01)      %   end of requested interval
        disp('==> End of time interval');
      end
      hold on;                                   % --- Accumulate trajectories
      for n=1:rays.NRAY
        ic=(n-1)*rays.odeDim;
        h=plot3(rays.YR(:,ic+rays.plt(1)),...
                rays.YR(:,ic+rays.plt(2)),...
                rays.YR(:,ic+rays.plt(3)));
        set(h,'LineWidth',sys.lineWidth);
      end
      hold off; view(0,0); 
      xlabel('major radius  R [m]');  zlabel('vertical direction  Z [m]');
      ylabel('radial wave vector k_R [m^{-1}]');

      if strcmp(rays.TYPE,'Amp')                      % --- Accumulate power
       for kray=1:rays.NRAY
        os=rays.odeDim*(kray-1); 
        nspec=size(plasma.amass,2);
        rr=rays.YR(:,os+1); zz=rays.YR(:,os+2);       % Ray coordinates in s
        rho=sqrt((rr-plasma.r0).^2 +zz.^2);
        tha=atan2(plasma.r0-rr,-zz)+pi/2;
        sray=solovev(rho,tha,plasma.r0,plasma.iaspr,plasma.elong,0);
        yr=[rays.YR(:,os+8) rays.YR(:,os+9+1:os+9+nspec)]; % log(E2) Pe Pi1 Pi2
        ns=size(plasma.depo,1); ds=1/ns;
        for ks=1:ns                              % Loop over radial intervals
         sl=(ks-1)*ds; sr=ks*ds;
         ind=find(sray>=sl & sray<=sr);
         if ~isempty(ind)
          for kk=1:length(ind)                   % Interpolate for bins in s
           im=max(ind(kk)-1,1 );          sm=sray(im);pm=yr(im,:);
           i0=    ind(kk)      ;          s0=sray(i0);p0=yr(i0,:);
           ip=min(ind(kk)+1,length(sray));sp=sray(ip);pp=yr(ip,:);
           if     (sm<sl)
            sray(i0)=sl; yr(i0,:)=yr(i0,:)+(p0-pm)*(s0-sl)./(s0-sm);
           elseif (sm>sr)
            sray(i0)=sr; yr(i0,:)=yr(i0,:)+(p0-pm)*(sr-s0)./(s0-sm);
           elseif (sp<sl)
            sray(i0)=sl; yr(i0,:)=yr(i0,:)+(pp-p0)*(sl-s0)./(sp-s0);
           elseif (sp>sr)
            sray(i0)=sr; yr(i0,:)=yr(i0,:)+(pp-p0)*(sr-s0)./(sp-s0);
           end
          end
          for kk=1:length(ind)-1                 % Accumulate in radial bins
           if ((ind(kk+1)-ind(kk))<2)
            plasma.depo(ks,:)=plasma.depo(ks,:) -yr(ind(kk),  2:end) ...
                                  +yr(ind(kk+1),2:end);
           end
          end %kk
         end %isempty
        end %ks
       end %kray
      end %if 'Amp'
  catch % this is something of a hack... all the recently calculated ray is lost
      err=lasterror;
      disp(err.message)
  end %try

%
%
% ----- Caustic: index of candidate rays -------------------------------------
elseif (strcmp(action,'caustic_which') && numel(rays.MON)>0)
  mon=sortrows(rays.MON,1);                      % Update history
  keep=6; jcaust=1;
  mon=mon(end-keep+1:end,:);                %   keep only most recent
  dt =diff(mon(:,1));
  dim1=size(mon,1);
  dim2=(size(mon,2)-1)/rays.NRAY;
  mon=reshape(mon(:,2:end)',dim2,rays.NRAY*dim1)';
  mon=reshape(mon(:,jcaust)',rays.NRAY,dim1)';   %   one col per ray
  dmon=diff(mon);
  dmondt=dmon./( dt*ones(1,size(dmon,2)) ).*rays.time./dim2;
  val=abs(mean(mon,1))';
  inc=abs(mean(dmondt,1))';
  %tune=0.05;
  tune=0.01;
  rays.caustic=find(abs(rays.monctc)<1E-8 | val<tune*inc);
  if (numel(rays.caustic)>0)                % Display result
    disp(['Time = ' num2str(rays.time)...
        '==> Signal caustic ray ' sprintf('#%i ',rays.caustic) ...
          ' -- Taylor |' sprintf('%0.3g| ',val) ...
                 '< f*|' sprintf('%0.3g| ',inc) ])
    disp(['                              Rays now in x=(' ...
                         sprintf(' %i',find(rays.inKspace==0)) ...
                 ') and k=(' sprintf(' %i',find(rays.inKspace==1)) ') space']);
  end
%
%
% ----- Caustic: proceed with Maslov transform -------------------------------
elseif (strcmp(action,'caustic_list') && ...
        numel(rays.MON)>0 && numel(rays.caustic)>0)
  dim=size(rays.y,1)/rays.NRAY;
  for kray=1:size(rays.caustic,2);          % Each ray independently
   k=rays.caustic(kray)-1;
   dWdr2=rays.y(k*dim+5);                   % Initial Theta derivatives
   dWdrz=rays.y(k*dim+6); 
   dWdz2=rays.y(k*dim+7);  
   GGTha=[dWdr2 dWdrz; dWdrz dWdz2];
   detGGTha=det(GGTha);
   invGGTha=inv(GGTha);
   rays.y(k*dim+5)=invGGTha(1,1);           % Transformed Theta derivatives
   rays.y(k*dim+6)=invGGTha(1,2);
   rays.y(k*dim+7)=invGGTha(2,2);
   idx=pi/4*sum(sign(eig(GGTha)));
   if (rays.inKspace(k+1))                  % Legendre transformation
    f=1./(2*pi*sqrt(abs(detGGTha)));sgn=+1; %   k -> x-space
   else                                    
    f=(2*pi)/sqrt(abs(detGGTha));   sgn=-1; %   x -> k-space
   end
   rays.y(k*dim+8)=rays.y(k*dim+8) + log(f^2);        %   amplitude
   rays.y(k*dim+9)=rays.y(k*dim+9) + idx + sgn* ...   %   phase
    (rays.y(k*dim+1)*rays.y(k*dim+3)+rays.y(k*dim+2)*rays.y(k*dim+4));
   rays.inKspace(k+1)=~rays.inKspace(k+1);
  end
  disp(['   ok, transformed ray ' sprintf('#%i ',rays.caustic) ...
        ' -- Rays now in x=(' sprintf(' %i',find(rays.inKspace==0)) ...
              ') and k=(' sprintf(' %i',find(rays.inKspace==1)) ') space']);
%
  if sum(rays.inKspace)                     % Short intervals in k-space
    timeintv=rays.timespan*0.05;
  else                                      % Larger in x-space
    timeintv=rays.timespan;
  end
  rays.timeintv=min(rays.timespan-rays.time,timeintv);
  rays.caustic=[];                          % Reset
%
%
% ----- Conversion: index of candidate rays ----------------------------------
elseif strcmp(action,'convert_which')

 if (numel(rays.MON)>0)
  mon=sortrows(rays.MON,1);                   % Sort the history by time.
  %mon = rays.MON;
  jconv=2;
  mon=mon(end-min(4,end-1):end,:);            % keep only most recent
  tim=mon(:,1);                               % time
  dim1=size(mon,1);dim2=(size(mon,2)-1)/rays.NRAY;
  mon=reshape(mon(:,2:end)',dim2,rays.NRAY*dim1)';
  TrD=reshape(mon(:,jconv)',rays.NRAY,dim1)';    % mon2 with one col per ray

  % rescale
  sca=tim(1);           
  tim=tim/sca; 
  %tt=rays.time/sca;

  A=[ones(size(tim)) tim tim.^2];           % quadratic fitting matrix

  lastwarn('');
  s=warning('off','MATLAB:rankDeficientMatrix');
  % cf=A\TrD;
  cf=A\abs(TrD);
  warning(s);

  [msgstr, msgid] = lastwarn;
  if strcmp(msgid,'MATLAB:rankDeficientMatrix')
      % the calculation of cf is suspect.  don't convert.
      rays.convert=[];
      disp('Problem with monitors during detection of conversion.');
  else
      %z0=cf(1,:)+cf(2,:)*tt+cf(3,:)*tt^2;
      %z1=2*tt*cf(3,:)+cf(2,:)
      z1=2*tim*cf(3,:)+cf(2,:);
      z2=2*cf(3,:);
      %tune=1;
      %convert=find(abs(z0)>tune*abs(z1) & z2<0);
      
      % added constraint that slope must change sign
      %rays.convert=find(abs(rays.moncnv)<1E-8 &...
      rays.convert=find(abs(rays.moncnv)<1E-8 &...
         z2'>0 &...
        (sign(z1(1))~= sign(z1(end))  ) ); % want min -> 2nd der >0 ?
  end

  if (numel(rays.convert)>0)                % Display result
    disp(['Time = ' num2str(rays.time)...
        '==> Signal convert ray ' sprintf('#%i ',rays.convert)]);
    %abs(rays.moncnv)
    %z2
    %2*tim*cf(3,:)+cf(2,:)
    % Seems to me (steve r) that if a conversion was detected, then the ray
    % should be updated so that it ends at z(t0) = z0.  Or this should be
    % done later (in convert_list)?  Or is it right already?
    
  end
  %abs(rays.moncnv)
  %z2
  %sign(z1(1))
  %sign(z1(end))
  %figure(3)
  %plot(tim,abs(TrD),'.-')%,tim,z1,'o-')
  %figure(1)
  %pause
  %diff(tim)
end
%
%
% ----- Conversion: proceed with mode-conversion -----------------------------
elseif (strcmp(action,'convert_list') && numel(rays.convert)>0)
  addRay=numel(rays.convert); rays.caustic=[];   % Initialize
  for kray=1:addRay                         % Each ray independently
   thisRay=rays.convert(kray);                   % Ray index
   k2x=rays.inKspace(thisRay)&strcmp(rays.TYPE,'Amp');% Temporarily back to x-space
   if (k2x)
     rays.caustic=thisRay; ray('caustic_list');
   end
   os=rays.odeDim*(thisRay-1);              % Calc velocity, acceleration
                                            %   and 2nd order derivatives
   z0 =rays.YR(end  ,os+1:os+rays.odeDim)'; time0=rays.TR(end);  
   zm1=rays.YR(end-1,os+1:os+rays.odeDim)'; tm1 = rays.TR(end-1);
   zm2=rays.YR(end-2,os+1:os+rays.odeDim)'; tm2 = rays.TR(end-2);
   dm2  =(tm2-tm1)*(tm2-time0); dm1=(tm1-tm2)*(tm1-time0); 
   d0=(time0-tm2)*(time0-tm1);
   zdot =zm2/dm2*(2*time0-time0-tm1) +zm1/dm1*(2*time0-time0-tm2)...
        +z0/d0*(2*time0-tm1-tm2);
   zddot=zm2/dm2* 2            +zm1/dm1* 2            +z0/d0* 2;
   
   mon =dispertok(0,z0,zdot,zddot,0,'Mon'); % Store yalf0
   zst0=dispertok(0,z0,zdot,zddot,0,'Sdl'); % Hyperbola saddle
   tol=1E-4; maxit=30;
   zst=zst0; zsto=zeros(size(zst)); it=0;
   while (norm((zst-zsto)./zst,inf)>tol && it<=maxit)
     zsto=zst; it=it+1;
     zst=dispertok(0,zst,zdot,zddot,0,'Sdl'); %zit=zst(1:4)  % outputs variable (debug)
   end
   dst=zst-zst0;
   if (it>maxit)
    disp('ray.m: hyperbola saddle point not converged--abort')%, pause
%   elseif (norm(dst(1:2))/plasma.r0>0.05 || norm(dst(3:4))/norm(plasma.kant)>2.)
    Estimate=zst0(1:4), Iterated=zst(1:4)   % outputs variable (debug)
    disp('ray.m: unlikely conversion point--abort')
   else
       
    % hack to save saddle point data...
    rays.zstList = [rays.zstList; zst];
       
    ytrs=dispertok(0,zst,zdot,zddot,0,'Trs'); % Transmitted ray
    
    rays.tau=[];  rays.beta=[];  % If addRay>1, this fails to save everything!
    ycnv=dispertok(0,zst,zdot,zddot,0,'Cnv'); % Converted ray: rays.tau and
                                              %   rays.beta set in this
                                              %   function call. 
    disp([' incm =' sprintf(' %0.3g',z0') ])
    disp([' conv =' sprintf(' %0.3g',ycnv') ])
    disp([' trsm =' sprintf(' %0.3g',ytrs') ])
    rays.y(os+1:os+rays.odeDim)=ycnv;                 % Update RHS vector
    rays.y=cat(1,rays.y,ytrs.'); rays.NRAY=rays.NRAY+1; %   append new ray
    
    % This line sets the space for the newly converted ray, but I've
    % changed how the code works, so do this elsewhere
    % rays.inKspace = [rays.inKspace 0];
    
    % is the new ray appended in the right dimension?  ie should it be
    % rays.y.' instead of ytrs.' ?

    if (k2x)                                % Switch back
     rays.caustic=thisRay; 
     ray('caustic_list');% Transform back
    end
   end
  end;

%
%
% ----- Conversion: matching and splitting -----------------------------------
elseif strcmp(action,'match_inc')
 global TM MM; TM=[]; MM=[];                % ----- Initialize
 keepInKspace=inKspace; T0=time;
 npts=31; scale=linspace(0,1,npts);
 adaptStepOptions = odeset('AbsTol',1E-6,'RelTol',1E-4,'Events','off');
 fixedStepOptions = odeset('AbsTol',1E-6,'RelTol',1E-4,'Events','off');
 for kk=2:npts
  tbck=[time time-0.10*scale(kk)*timeintv]; % 1st backward from sigma_0
  tfwd=[time time+0.05*scale(kk)*timeintv]; % 1st forward  from sigma_0
%  tbck=[time time-0.08*scale(kk)*timeintv]; % 2nd backward from sigma_0
%  tfwd=[time time+0.10*scale(kk)*timeintv]; % 2nd forward  from sigma_0
%
  addRay=prod(size(convert));
  for kray=1:addRay                         % Each ray independently
   inKspace=keepInKspace(kray);
%
   thisRay=convert(kray);                   % ----- Incident+transmitted ray
   os=odeDim*(thisRay-1);
   zi=y(os+1:os+odeDim);                    % Backwards starting from z0
   [ti,zi]=feval('ode23s','match',tbck,zi,adaptStepOptions);
   dt=0.5*(ti(end-1)-ti(end)); zi=zi(end,:);
   tdiff=ti(end):dt:ti(end-1);              % Velocity and Acceleration
   [ti,zi]=feval('ode23s','match',tdiff,zi,fixedStepOptions);
   zidiff=diff(zi,1); zi=zi(2,:); ti=ti(2);
   zidot =0.5*sum(zidiff,1)/dt;
   ziddot=diff(zidiff,1)/(dt*dt);

   mon =dispertok(0,zi,zidot,ziddot,0,'Mon');% Store matching position
   zst0=dispertok(0,zi,zidot,ziddot,0,'Sdl');% Hyperbola saddle
   tol=1E-4; maxit=30;
   zst=zst0; zsto=zeros(size(zst)); it=0;
   while (norm((zst-zsto)./zst,inf)>tol & it<=maxit)
     zsto=zst; it=it+1;
     zst=dispertok(0,zst,zidot,ziddot,0,'Sdl');
     zit=zst(1:4);
   end
   dst=zst-zst0;
   if (it>maxit)                            % Check if saddle converged
    disp(['ray.m: saddle point not converged']), pause
   elseif (norm(dst(1:2))/r0>0.05 || norm(dst(3:4))/norm(kant)>2.)
    Estimate=zst0(1:4), Iterated=zst(1:4);
    disp(['ray.m: unlikely conversion point']), pause
   else
    zt=dispertok(0,zi, zidot,ziddot,0,'Mon')';% Store inc position
    m1=dispertok(0,zi, zidot,ziddot,0,'Mch')';%   eta,tau from inc pos
    m2=dispertok(0,zst,zidot,ziddot,0,'Mch')';%   eta,tau from z*
    zt=dispertok(0,zst,zidot,ziddot,0,'Trs')';% Transmitted ray
    zc=dispertok(0,zst,zidot,ziddot,0,'Cnv')';% Converted ray
    TM=[TM; ti];
    MM=[MM; [zit(1) zi(1) zt(1) sqrt(exp(zi(8))) sqrt(exp(zt(8))) m1 m2]];
   end
                                            % ----- Converted ray
   zi=y(os+1:os+odeDim);                    % Forward from z0
   [tc,zc]=feval('ode23s','match',tfwd,zi,adaptStepOptions);
%
   dt=0.5*(tc(end)-tc(end-1)); zc=zc(end,:);
   tdiff=tc(end-1):dt:tc(end);              % Velocity and Acceleration
   [tc,zc]=feval('ode23s','match',tdiff,zc,fixedStepOptions);
   zcdiff=diff(zc,1); zc=zc(2,:); tc=tc(2);
   zcdot =0.5*sum(zcdiff,1)/dt;
   zcddot=diff(zcdiff,1)/(dt*dt);
%
   mon =dispertok(0,zc,zcdot,zcddot,0,'Mon');% Store matching position
   zst0=dispertok(0,zc,zcdot,zcddot,0,'Sdl');% Hyperbola saddle
   tol=1E-4; maxit=30;
   zst=zst0; zsto=zeros(size(zst)); it=0;
   while (norm((zst-zsto)./zst,inf)>tol & it<=maxit)
     zsto=zst; it=it+1;
     zst=dispertok(0,zst,zcdot,zcddot,0,'Sdl');
     zit=zst(1:4);
   end
   dst=zst-zst0;
   if (it>maxit)                            % Check if saddle converged
    disp(['ray.m: saddle point not converged']), pause
   elseif (norm(dst(1:2))/r0>0.05 || norm(dst(3:4))/norm(kant)>2.)
    Estimate=zst0(1:4), Iterated=zst(1:4);
    disp(['ray.m: unlikely conversion point']), pause
   else
    zt=dispertok(0,zc, zcdot,zcddot,0,'Mon')';% Store inc position
    m1=dispertok(0,zc, zcdot,zcddot,0,'Mch')';%   eta,tau from inc pos
    m2=dispertok(0,zst,zcdot,zcddot,0,'Mch')';%   eta,tau from z*
    zt=dispertok(0,zst,zidot,ziddot,0,'Trs')';% Transmitted ray
    zc=dispertok(0,zst,zcdot,zcddot,0,'Cnv')';% Converted ray
    TM=[tc; TM];
    MM=[[zit(1) zc(1) zt(1) sqrt(exp(zc(8))) sqrt(exp(zt(8))) m1 m2];MM];
   end
  end
 end
 mid=size(MM,1)/2+1;
 TM=TM*1E6; MM0=MM(mid,:);
 MM0=ones(prod(size(TM)),1)*MM0;
 MMR=(MM-MM0)./MM0, figure(4);
%
 Match=cat(2,TM,MM)                         % Output results
 h=plot(TM,MMR(:,2),'b',TM,MMR(:,4)/10,'b--', ...
        TM,MMR(:,3),'r',TM,MMR(:,5)/10,'r--',TM,MMR(:,6),'g',TM,MMR(:,7),'k');
 set(h,'LineWidth',lineWidth); hold on
 legend('R_{1m} or R''_{1m}','|E_{1m}|/10 or |E''_{1m}|/10', ...
        'R_{2m}','|E_{2m}|/10', ...
        '\tau(z_{1m})','\tau(z_*)', ...
        'Location','SouthEast')
%%%%%
 lbl=['Relative deviation from values matched at t_0=', ...
      sprintf('%4.3f',T0*1E6)];
 ylabel(lbl)
 xlabel('Matching time t_{1m}<t_0 and t''_{1m}>t_0 [\mu sec]')
 axis([TM(end) TM(1) -0.025 0.025]); hold off;
 inKspace=keepInKspace;
%
%
% ----- Accumulate output in plot --------------------------------------------
elseif strcmp(action,'accumul')
  hold on;                                 % accumulates runs
  for n=1:size(rays.YR,2)/rays.odeDim
    jc=(n-1)*rays.odeDim; cl=sys.pltCol(mod(n-1,size(sys.pltCol,2))+1);
    h=plot3(rays.YR(:,jc+rays.plt(1)),...
            rays.YR(:,jc+rays.plt(2)),...
            rays.YR(:,jc+rays.plt(3)),cl);
    set(h,'LineWidth',sys.lineWidth);
  end
  hold off; view(0,0);
  xlabel('major radius  R [m]');  
  zlabel('vertical direction  Z [m]');
  if (rays.odeDim==4)
    if     (rays.plt(2)==3) ylabel('radial wave vector k_R [m^{-1}]');
    elseif (rays.plt(2)==4) ylabel('vertical wave vector k_Z [m^{-1}]');
    end;
  elseif (rays.odeDim==6)
    if     (rays.plt(2)==4) ylabel('radial wave vector k_R [m^{-1}]');
    elseif (rays.plt(2)==6) ylabel('vertical wave vector k_Z [m^{-1}]');
  end;
  end
%
%
% ----- History: monitors ----------------------------------------------------
elseif (strcmp(action,'history1') && ~strcmp(rays.TYPE,'Trj'))
  rays.MON=sortrows(rays.MON,1);                      % Update history
  dim1= size(rays.MON,1);
  if isempty(rays.MONNorm)                       % Remember initial normalization
    rays.MONNorm=mean(rays.MON,1);
  end
  tmp=rays.MON(:,2:end);                         % Plot curves
  tim=reshape(ones(rays.NRAY,1)*rays.MON(:,1)',rays.stp*rays.NRAY,1);
  nbr=reshape((1:rays.NRAY)'*ones(1,rays.stp) ,rays.stp*rays.NRAY,1);
  dat=reshape(tmp',size(tmp,2)/rays.NRAY,size(tmp,1)*rays.NRAY)';
  srt=sortrows(cat(2,tim,nbr,dat),[2 1]);   %   historical order
  srt(:,1)=srt(:,1)*1E6;                    %    1 time in micro-seconds
                                            %    2 ray number
  srt(:,3)=srt(:,3)./rays.MONNorm(2);       %    3 norm caustics monitor
  srt(:,4)=srt(:,4)./rays.MONNorm(3);       %    4      conversion monitor
                                            %    5 coupling coef. eta^2
  if (size(tim,1)<10)
    RAY_monitors_t_nbr_m1_m2_exp= ...
      cat(2,srt(:,1)./max(srt(:,1)),...
            srt(:,2:4),exp(-pi*srt(:,5)))
  end;
  t=1; nr=2; m1=3; m2=4; e2=5; Rc=6; kr=7;
  plot(srt(:,t),srt(:,m1),'k.',...
       srt(:,t),srt(:,m2),'r.',...
       srt(:,t),exp(-pi*srt(:,e2)),'b.' )
  axis([srt(1,t) srt(end,t) -0.5 1.5])
  xlabel('time [\mu sec]'); ylabel('monitored quantities [a.u.]');
  legend('4|\nabla\theta|-det(\nabla\nabla\theta)^{-1/2}',...
         'Tr(D)','exp(-\pi\eta^2)',0);
%
%
% ----- History: conversion parameters ---------------------------------------
elseif (strcmp(action,'history2') && ~strcmp(rays.TYPE,'Trj'))
  tim =rays.MON(:,1)*1E6;   cstc=rays.MON(:,2);    conv=rays.MON(:,3);
  eta2=rays.MON(:,4);       Rc  =rays.MON(:,5);    kr  =rays.MON(:,6);
  plot(tim,Rc, tim,kr, tim,exp(-pi*eta2),'.');
  xlabel('time [\mu sec]'); ylabel('monitored quantities');
  legend('R_c','kr_c','\tau');
%
%
% ----- History: ODE variables -----------------------------------------------
elseif strcmp(action,'history3')
  if (rays.odeDim==4)                            % Nbr subplots
   vp=2;
  else
   vp=4;
  end
  dim=size(rays.YR,2)/rays.NRAY;
  for kk=0:rays.NRAY-1
   off=kk*dim; cl=sys.pltCol(mod(kk,size(sys.pltCol,2))+1);
   subplot(vp,2,1);hold on;h=plot(rays.TR*1E6,rays.YR(:,off+1)-plasma.r0,cl);
                          set(h,'LineWidth',sys.lineWidth); ylabel('R-R0 [m]')
   subplot(vp,2,2);hold on;h=plot(rays.TR*1E6,rays.YR(:,off+2),cl);
                          set(h,'LineWidth',sys.lineWidth); ylabel('Z [m]')
   subplot(vp,2,3);hold on;h=plot(rays.TR*1E6,rays.YR(:,off+3),cl);
                          set(h,'LineWidth',sys.lineWidth); ylabel('k_R [1/m]')
   subplot(vp,2,4);hold on;h=plot(rays.TR*1E6,rays.YR(:,off+4),cl);
                          set(h,'LineWidth',sys.lineWidth); ylabel('k_Z [1/m]')
   if (rays.odeDim==4)
     xlabel('time [\mu sec]'); subplot(vp,2,3); xlabel('time [\mu sec]');
   else
     %supress=(1-rays.inKspace(kk+1)+1E-6); 
     supress = 1;
     %nK=sqrt(rays.YR(:,off+3).^2+plasma.kant(2)^2+rays.YR(:,off+4).^2); % |k|^{-1}
     nK=ones(size(supress));
     subplot(4,2,5);hold on;h=plot(rays.TR*1E6,supress*rays.YR(:,off+5)./nK,cl);
             set(h,'LineWidth',sys.lineWidth); ylabel('\partial_{RR}\Theta [1/m^2]')
     subplot(4,2,6);hold on;h=plot(rays.TR*1E6,supress*rays.YR(:,off+7)./nK,cl);
             set(h,'LineWidth',sys.lineWidth); ylabel('\partial_{ZZ}\Theta [1/m^2]')
     subplot(4,2,7);hold on;h=plot(rays.TR*1E6,supress*sqrt(exp(rays.YR(:,off+8))),cl);
             set(h,'LineWidth',sys.lineWidth); ylabel('|E|'); 
             xlabel('time [\mu sec]');
     subplot(4,2,8);hold on;h=plot(rays.TR*1E6,supress*rays.YR(:,off+9)/(2*pi),cl); 
             set(h,'LineWidth',sys.lineWidth); ylabel('\Theta [rad/2\pi]');
             xlabel('time [\mu sec]');
   end
  end
%
%
% ----- History: power deposition --------------------------------------------
elseif strcmp(action,'history4')
  ns=size(plasma.depo,1); ds=1/ns;
  subplot(2,1,1); hold on
  dim=size(rays.YR,2)/rays.NRAY;
  for kk=0:rays.NRAY-1
    off=kk*dim; cl=sys.pltCol(mod(kk,size(sys.pltCol,2))+1);
    h=plot(rays.TR*1E6,rays.YR(:,off+10),[cl '--'], ...
           rays.TR*1E6,rays.YR(:,off+11),[cl '-.'], ...
           rays.TR*1E6,rays.YR(:,off+12),[cl ':']);
    set(h,'LineWidth',sys.lineWidth);
  end
  ylabel('power absorption  W(t)')
  xlabel('time [\mu sec]');

  subplot(2,1,2); cla reset; s=ds:ds:1; s(1)=0;
  h=plot(s,plasma.depo(:,1),'k--',s,plasma.depo(:,2),'k-.',s,plasma.depo(:,3),'k:');
  set(h,'LineWidth',sys.lineWidth); ylabel('power deposition  P(s)')
  axis([0 1 0 max(max(plasma.depo))])

% prof=[zeros(size(amass));depo];
% for kk=3:ns+1
%  prof(kk,:)=prof(kk,:)+prof(kk-1,:);
% end
% tot=sum(prof(end,:)); prof=prof/tot;
% subplot(2,1,2); s=0:ds:1;
%    h=plot(s,prof(:,1),'b--',s,prof(:,2),'r-.',s,prof(:,3),'g:');
%    set(h,'LineWidth',lineWidth); ylabel('P_{tot}^{-1} \int_0^s P(r) dr')
  xlabel('normalized radius  s=(\psi/\psi_s)^{1/2}')
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% GUI -- Graphical User Interface events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% ----- GUI event: change data set -------------------------------------------
elseif strcmp(action,'newData')
  figNumber = watchon;   h = get(figNumber,'Children');
  dataHndl= findobj(h,'flat','Tag','data');
  s=get(dataHndl,'Value');
  type=dataTypes(s,:);
  cmdStr=data(type);                        % Read parameters from data.m
  set(windowHndl(2),'String',cmdStr);       % Activate scenario
  ray('evalText');                          % Input user data from textfield
%
%
% ----- GUI event: editable text field ---------------------------------------
elseif strcmp(action,'evalText')
  evalmcw(windowHndl(2));                   % evaluate content
  ray('auxval');                            % Auxiliary data
  ray('start');                             % New Ray IC
  ray('newWindow');                         % Compute new plot
%
%
% ----- GUI event: update plot window ----------------------------------------
elseif strcmp(action,'newWindow')
  figNumber = watchon;  h = get(figNumber,'Children');
  hndl= findobj(h,'flat','Tag','window');
  s=get(hndl,'Value');
  type=windowTypes(s,:);
%
  if strcmp(type,'Parameters  ');           % Toggle windows
   set(windowHndl(1),'Visible','off');
   set(windowHndl(2),'Visible','on');
   axis off;
  else
   set(windowHndl(2),'Visible','off');
   set(windowHndl(1),'Visible','on');
   axis on;
  end
%
  if strcmp(type,'Parameters  ');           % Actions depending on selector
    ray('auxval');
  elseif strcmp(type,'EquilData 1D');
    cla reset; ray('equil1D');
  elseif strcmp(type,'Roots 1D    ');
    ray('roots');
  elseif strcmp(type,'Mesh 2D     ');
    ray('mesh'); 
  elseif strcmp(type,'EquilData 2D');
    cla reset; ray('equil2D');
  elseif strcmp(type,'Roots 3D    ');
    cla reset; ray('roots3D');
  elseif strcmp(type,'RayInitializ');
    cla reset; ray('start');
  elseif strcmp(type,'RayPropagate');
    ray('propagate')
  elseif strcmp(type,'Hist monitor');
    legend off; cla reset; ray('history1')
  elseif strcmp(type,'Hist ODE var');
    legend off; cla reset; ray('history2')
  end
%
%
% ----- GUI event: print plot to file ----------------------------------------
elseif strcmp(action,'print')
 fn=cat(2,'Fig',datestr(now,13),'.tif')
 print('-dtiff','-r80',fn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Error reports %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
else
  disp(['Error: cannot perform ray(''',action,''')'])

  if (strcmp(action(1:4),'hist') & strcmp(TYPE,'Trj'))
    disp(['       set TYPE=''Con'' or ''Amp'' for conversion monitoring'])
  end

  if (strcmp(action(1:4),'conv') & strcmp(TYPE,'Trj'))
    disp(['       set TYPE=''Con'' or ''Amp'' for ray conversion'])
  end

  if (strcmp(action(1:4),'caus') & ~strcmp(TYPE,'Amp'))
    disp(['       set TYPE=''Amp'' for Legendre transformation'])
  end

  if (strcmp(action(1:4),'caus') & prod(size(caustic))==0)
    disp(['       set caustic=[k l m] to specify which rays to transform'])
  end

  if (strcmp(action(1:4),'caus') & prod(size(MON))==0)
    disp(['       First need to propagate a few steps...'])
  end

  if (strcmp(action(1:4),'conv') & prod(size(convert))==0)
    disp(['       set convert=[k l m] to specify which rays to convert'])
  end
end;

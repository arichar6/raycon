function [tspan,y,options] = inittok(tintv,y)
% ===== Initial conditions ===================================================
global plasma rays cnst

% ----- ODE integration parameters
if strcmp(rays.TYPE,'Trj')         
  rays.evts='off'; rays.plt=[1 3 2]; rays.odeDim=4;
  %atol = 1E-10*ones(1,rays.odeDim);
  atol = [1E-8 1E-8 1E-6 1E-6];
elseif strcmp(rays.TYPE,'Con')
  rays.evts='on';  rays.plt=[1 3 2]; rays.odeDim=4;
  atol = 1E-7*ones(1,rays.odeDim);
elseif strcmp(rays.TYPE,'Amp')
  rays.evts='on';  rays.plt=[1 3 2]; rays.odeDim=9+length(plasma.amass);
  % Specify vector for tolerance, so ray position is more accurate than
  % other quantities.
  atol = [1E-8*ones(1,4) 1E-2*ones(1,rays.odeDim-4)];
else
  error(['trajectory: Unknown TYPE ''' rays.TYPE '''.']);
end;
%options = odeset('AbsTol',1E-7,...
%                 'RelTol',1E-6,...
%options = odeset('AbsTol',1E-6,...
%                 'RelTol',1E-4,...
%options = odeset('AbsTol',atol,... 
%                 'RelTol',1E-6,...
options = odeset('AbsTol',atol,...
                 'RelTol',1E-6,...
                 'Events',rays.evts,...
                 'InitialStep',rays.initialstep,...
                 'OutputFcn','odephas3',...
                 'Vectorized','off',...
                 'OutputSel',rays.plt,...
                 'Refine',16);
%
if (numel(y)>0)                       % ----- Resume from previous data
    % Reset the ray so it is on the dispersion surface?
%     y0=y;                                % initial guess
%     if (rays.odeDim==4)                  % Wave vector
%       kdir=[0 0 1 0]';                   %   solve kr=kr(kz)
%     elseif (rays.odeDim>=9)              % Wave vector, amplitude & phase
%       kdir=[0 0 1 0 0 0 0 0 0 zeros(size(plasma.amass))]';
%     end;
%     for i=1:rays.NRAY % Adjust to start on the dispersion surface
%         %norm(y0(2:3,i))
%         kcor= fzero('dispertok',y0(3,i), ...
%                     optimset('Display','off'), ...
%                     y0(:,i),kdir,zeros(size(kdir)),0,'Dsp');
%         y0(3,i)=y0(3,i)+kcor;
%         rays.kray0(i,1)=rays.kray0(i,1)+kcor;
%         kcor= fzero('dispertok',y0(3,i), ...    
%                     optimset('Display','off'), ...
%                     y0(:,i),kdir,zeros(size(kdir)),0,'Dsp');
%         y0(3,i)=y0(3,i)+kcor;
%         if (rays.odeDim>=9)                         % Focusing
%           y0(5:9,i)=dispertok(0,y0(:,i),y0(:,i),y0(:,i),0,'Ant');
%         end
%     end
%     y = reshape(y0,rays.odeDim*rays.NRAY,1);
else                                  % ----- Set Initial conditions
  [rho,r,z]=mapFlux(rays.sray0,rays.thray0);         % Positions in (R,Z)
  zro=zeros(size(r')); %one=ones(size(r'));
  rays.kray0=ones(rays.NRAY,1)*plasma.kant;
  if (rays.odeDim==4)                                % Wave vector
    y0=cat(2,r',z',rays.kray0(:,1),rays.kray0(:,3))';%   initial guess
    kdir=[0 0 1 0]';                                 %   solve kr=kr(kz)
  elseif (rays.odeDim>=9)                            % Wave vector, amplitude & phase
    y0=cat(2,r',z',rays.kray0(:,1),rays.kray0(:,3),...
             zro,zro,zro,zro,zro,...
             zeros(rays.NRAY,size(plasma.amass,2)))';
    kdir=[0 0 1 0 0 0 0 0 0 zeros(size(plasma.amass))]';
  end;
%end
  
  for i=1:rays.NRAY                             % Adjust to start on MS surface
    T=dispertok(0,y0(:,i),y0(:,i),y0(:,i),0,'Msw');
    kms=sqrt((plasma.omega/cnst.c)^2*T);        % Estimate for start of iteration
    kcor= fzero('dispertok',kms, ...
                optimset('Display','off'), ...
                y0(:,i),kdir,zeros(size(kdir)),0,'Dsp');
    y0(3,i)=y0(3,i)+kcor;
    rays.kray0(i,1)=rays.kray0(i,1)+kcor;
    kms=sqrt((plasma.omega/cnst.c)^2*T);% iterate the root finding twice,
    kcor= fzero('dispertok',kms, ...    %   but I'm not sure why.
                optimset('Display','off'), ...
                y0(:,i),kdir,zeros(size(kdir)),0,'Dsp');
    y0(3,i)=y0(3,i)+kcor;
    rays.kray0(i,1)=rays.kray0(i,1)+kcor;
    if (rays.odeDim>=9)                         % Focusing
      y0(5:9,i)=dispertok(0,y0(:,i),y0(:,i),y0(:,i),0,'Ant');
    end
  end;
  y = reshape(y0,rays.odeDim*rays.NRAY,1);
end; % if numel(y) > 0

tspan=[rays.time rays.time+tintv];                        % ----- New time interval
rays.ht=zeros(1,6);           rays.ht(1)=rays.time;       % reset 6 levels history
rays.hy=zeros(rays.odeDim*rays.NRAY,6); rays.hy(:,1)=y;   %   positions

if (strcmp(rays.TYPE,'Con')||strcmp(rays.TYPE,'Amp'))     % Eikonal phase correction
  rays.pol=dispertok(0,y,y,y,0,'Pol');                    %   polarizations 
  rays.polDim=numel(rays.pol);
  rays.hp=zeros(rays.polDim,6); rays.hp(:,1)=reshape(rays.pol,rays.polDim,1);
  rays.poldot=zeros(rays.polDim,rays.NRAY);               %   derivative
else
  rays.poldot=0;
end
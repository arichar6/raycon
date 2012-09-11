%% RayCon Cells 
%  A set of MatLab cells which allow interaction with the RayCon
%  calculations.  Designed as an interactive alternative to the sequence
%  programming in ray.m.

%% Initialize global variables
global ...     % Plasma configuration
    plasma...  % Struct containing the physical configuration variables
    rays...    % Strcut with data about rays
    cnst...    % Struct with physical constants
    sys        % Struct with system preferences (lineWidth, etc)

plasma = initPlasma;
rays   = initRays;
cnst   = initCnst;
sys    = initSys;

%% Set system parameters for CMod ICRF calculation

data('cmod');                              % Load the CMod paramters
rays.TYPE='Con';                           % type of tracing to do
%rays.TYPE='Trj';                           % type of tracing to do
rays.odeDim = 4; % z and S matrix
rays.timespan = 5e-2;                       % Time span of evolution [1/freq]
rays.aThres   = 0.05;                      % Amplitude threshold
rays.NRAY = 1;                             % number of rays to trace
rays.time = 0;

plasma.thant = [-.5;.5];
th_oneray=0.001;
plasma.kant = [-31.5 -10 0.];

%% Set initial conditions for some rays
% rays.sray0,thray0; plasma.kant determine ic's in ray('start')/inittok
% numerical tolerance for ode solver also set in inittok

numRays = 1;  % number of rays to start (since I don't know how NRAY should work)
rays.RayICList = [];
rays.RayIniTimeList = [];
rays.RayTRList = [];
rays.RayYRList = [];
rays.RayIndList = [];
rays.RayConvList = [];
rays.RayEndList = [];
rays.RayInKList = [];
rays.InK = [];
rays.y = [];
if numRays == 1
    tmp = plasma.thant(1):...
       (plasma.thant(2)-plasma.thant(1))/(19):plasma.thant(2);
    %th0_list = tmp(11);
    th0_list = th_oneray;
else
    th0_list = plasma.thant(1):...
       (plasma.thant(2)-plasma.thant(1))/(numRays-1):plasma.thant(2);
end
for th0=th0_list
    rays.sray0 = plasma.sant;
    rays.thray0 = th0;
    [rho,r,z]=mapFlux(rays.sray0,rays.thray0);         % Positions in (R,Z)
    tspan=[rays.time rays.time+rays.timespan];
    rays.tspan = tspan;
    % ODE options are set in inittok
             
    % set IC's
    y = adjust_disp_m([r,z,plasma.kant(1:2:3)],0);  %m=last argument
%     if (rays.odeDim>=9)                             % Focusing
%         y0=[y.' zeros(1,rays.odeDim - 4)];
%         y0=dispertok(0,y0,y0,y0,0,'Ant');
%         y=[y.' y0.' zeros(1,rays.odeDim - 9)].';
%     end
    rays.y=y;
    rays.RayICList = [rays.RayICList; y.'];
    rays.RayIniTimeList = [rays.RayIniTimeList; 0];
        % each row is a new ray to trace
end


%% Plot zeros of dispersion surface in the mid-plane
if isempty(who('R'))
    [r_small,z_small] = meshgrid(0.405:.005:0.9,0.001);
    k=linspace(-2*abs(plasma.kant(1)),2*abs(plasma.kant(1)),100);
    R=zeros(size(r_small,1),size(r_small,2),size(k,2)); 
    Z=R; K=R; D=R; % Vec_R=R; Vec_kR=R;
    r_plane =reshape(r_small,numel(r_small),1);
    z_plane =reshape(z_small,numel(z_small),1);

    tmpDim = rays.odeDim;
    rays.odeDim = 4;
    for i3=1:size(k,2)
        R(:,:,i3)=r_small; Z(:,:,i3)=z_small; 
        K(:,:,i3)=ones(size(r_small)).*k(i3);
        tmp =cat(2,r_plane,z_plane,ones(size(r_plane))*k(i3),...
           ones(size(r_plane))*plasma.kant(3));
        yv=reshape(tmp',1,numel(tmp))';
        %Dv=dispertok(0.,yv,yv,yv,0,'Eig');
        %Dv=dispersion_new(yv,'Dsp');
        for i4=1:4:length(yv)
            Dv=disp_eig(yv(i4:i4+3),'Dsp');
            %D(:,i4,i3)=reshape(Dv',size(r_small,1),size(r_small,2));
            D(:,1+(i4-1)/4,i3)=Dv;
%             Dv=disp_eig(yv(i4:i4+3),'Trj');
    %         tmp=reshape(Dv',size(r_small,1),size(r_small,2),4);
    %         Vec_R(:,:,i3)=tmp(:,:,1);
    %         Vec_kR(:,:,i3)=tmp(:,:,3);
%             Vec_R(:,1+(i4-1)/4,i3)=Dv(1:4:end);
%             Vec_kR(:,1+(i4-1)/4,i3)=Dv(3:4:end);
        end
    end
    rays.odeDim = tmpDim;
end

%%
figure(2)
contour(squeeze(R),squeeze(K),squeeze(D),[0 0],'k');

%%
% figure(2)
% %surf(squeeze(R),squeeze(K),squeeze(Vec_R./(abs(Vec_R)<100)))
% %surf(squeeze(R),squeeze(K),squeeze(Vec_kR./(abs(Vec_kR)<1e6)))
% lim = (abs(Vec_R)<1e7) & (abs(Vec_kR)<1e9);
% v1=squeeze(Vec_R./lim);
% v2=squeeze(Vec_kR./lim);
% v1(~isfinite(v1))=0;
% v2(~isfinite(v2))=0;
% hold on
% quiver(squeeze(R),squeeze(K),v1,v2,'ShowArrowHead','off')
% hold off

%% Now it's time to try some rays!

% turn on debug output
%sys.debugMon = 1;
figure(2)
contour(squeeze(R),squeeze(K),squeeze(D),[0 0],'k');
j = 1;
while j<=size(rays.RayICList,1) % loop through rays
    beep % alert user that new ray is starting...
    % Set the position part of IC's for the next ray
    rays.y = rays.RayICList(j,:).';
    rays.time = rays.RayIniTimeList(j);
    rays.inKspace = 0;
    % Reset monitored quantities
    numConvert=0; % limit ray tracing time by number of conversion... 
                  % ...rough estimate of lost power
    rays.MON=[]; rays.MONNorm=[]; rays.TR=[]; rays.YR=[];
    rays.monErrAbort = 0;
    rays.stp=0; 
    rays.inKspace=zeros(1,rays.NRAY);
    rays.timeintv=rays.timespan;
    %rays.initialstep = 1e-4*rays.timespan;
    rays.initialstep = 1e-7*rays.timespan;
    disp(['Ray#' sprintf('%i  init =',j)...
                 sprintf(' %0.3g',rays.y(1:rays.odeDim)) ])
    figure(2)
    hold on
    plot(rays.y(1),rays.y(3),'o','lineWidth',1)
    hold off
    title(['Time = ' num2str(rays.time)])
    
    while (rays.time < rays.timespan && ~rays.monErrAbort ...
            && numConvert<3)
        if rays.time ~=0 && length(rays.TR)>2
            rays.initialstep = abs(rays.TR(end) - rays.TR(end-2));
        end
        if rays.initialstep == 0 % if the step is too small,  
            break;               %  there was some error in tracing, so stop.
        end
        ray('propagate')
        rays.RayInKList = [rays.RayInKList; rays.InK];
        figure(2)
        hold on
        plot(rays.YR(:,1),rays.YR(:,3),'-','lineWidth',2)
        hold off
        title(['Time = ' num2str(rays.time)])
        %rays.convert
        %rays.caustic
        %pause
        
        % Deal with mode conversion
        if (numel(rays.convert)~=0)
            ray('convert_list')
            % save ic for transmitted ray and data about conversion
            if size(rays.y,2)>1 % the conversion calc worked, save it
                numConvert = numConvert+1; % update conversion counter
                % ConvList = [time, inc ray num, 
                %     trans ray num, tau, beta, saddle point]
                rays.RayConvList = [rays.RayConvList;  rays.time, j, ...
                  size(rays.RayICList,1)+1, rays.tau, rays.beta, rays.y(1,:)];
                rays.RayICList(end+1,:) = rays.y(2,:).';
                rays.RayIniTimeList(end+1) = rays.time;
                hold on
                plot(rays.y(1,1),rays.y(1,3),'rx')              % saddle point
                plot(rays.y(2,1),rays.y(2,3),'bx')              % transmitted ray
                plot(rays.YR(end,1),rays.YR(end,3),'kx')        % converted ray
                title(['Time = ' num2str(rays.time)])
                hold off
            end
            rays.y = rays.YR(end,:).';
            rays.NRAY = 1;
            rays.convert = [];
            % rays.MON=[]; rays.TR=[]; rays.YR=[]; rays.stp=0;
        end
        % Deal with caustic
        if numel(rays.caustic)~=0
%             figure(3)
%             ray('history3')
             reply = input('Switch to/from k-space? Y/N [Y]: ', 's');
             if isempty(reply)
                 reply = 'Y';
             end
             if reply == 'Y'
                ray('caustic_list')
             else
                 disp('Switching canceled.');
             end
            % pause
        end
    end

    % Save the traced data for this ray
    % In some cases, the tracing fails, and the data is lost -> TR and YR
    % will be [].  The following would then fail...
%     if ~isempty(rays.TR) % try to recover if the tracing failed
%         rays.YR=rays.y.';
%         rays.TR=rays.time;
%         rays.monErrAbort=1; % make sure the loop stops
%     end
    rays.RayIndList = [rays.RayIndList; size(rays.RayTRList,1)+1];
    rays.RayTRList  = [rays.RayTRList; rays.TR];
    rays.RayYRList  = [rays.RayYRList; rays.YR];
    j = j+1;
    
%     figure(3)
%     clf
%     ray('history3')
    %pause
end  

rays.RayEndList = [rays.RayIndList(2:end)-1; length(rays.RayTRList)];


%% Plot all the rays together
figure(2)
clf
contour(squeeze(R),squeeze(K),squeeze(D),[0 0],'k');
hold on
% plot the rays
for j = 1:size(rays.RayIndList,1)
%for j = 1:3
    plot(rays.RayYRList(rays.RayIndList(j):rays.RayEndList(j),1),...
         rays.RayYRList(rays.RayIndList(j):rays.RayEndList(j),3),...
         'b-','lineWidth',2)
    % plot ic of the rays
    plot(rays.RayYRList(rays.RayIndList(j),1),...
         rays.RayYRList(rays.RayIndList(j),3),'bo')
end

%% 
if isempty(who('D2'))
%    r_line=0.2:.001:1;
%    z_line=-.5:.005:.5;
    r_line=0.4:.005:.9;
    z_line=-.4:.01:.4;
    [r_small,z_small] = meshgrid(r_line,z_line);
    R2=zeros(size(r_small,1),size(r_small,2),1); 
    Z2=R2; K2=R2; D2=R2; Vec_R2 = R2; Vec_Z2 = R2;
    
    k0=rays.RayICList(3:4);
    tmpDim = rays.odeDim;
    rays.odeDim = 4;
    for i1=1:length(r_line)
            tmp=[r_line(i1)+zeros(size(z_line))',z_line',...
               zeros(size(z_line))',zeros(size(z_line))'];
            yv=reshape(tmp',1,numel(tmp))';
            %Dv=dispertok(0.,yv,yv,yv,0,'Dsp');
            Dv=disp_eig(yv,'Dsp');
            D2(:,i1)=Dv;

%             tmp=[r_line(i1)+zeros(size(z_line))',z_line',...
%                k0(1)*ones(size(z_line))',k0(2)*ones(size(z_line))'];
%             yv=reshape(tmp',1,numel(tmp))';
%             Dv=disp_eig(yv,'Trj');
%             Vec_R2(:,i1)=Dv(1:4:end);
%             Vec_Z2(:,i1)=Dv(2:4:end);
    end
    rays.odeDim = tmpDim;
end

%% Plot in real space
figure(2)
clf 
plasma=plotFlux(plasma);
hold on
contour(r_line,z_line,squeeze(D2(:,:,1)),[0 0],'k');
%%
% lim = (abs(Vec_R2)<5e3) & (abs(Vec_Z2)<5e3);
% v21=squeeze(Vec_R2./lim);
% v22=squeeze(Vec_Z2./lim);
% v21(~isfinite(v21))=0;
% v22(~isfinite(v22))=0;
% quiver(r_small,z_small,v21,v22,'ShowArrowHead','off')

% plot the rays
for j = 1:size(rays.RayIndList,1)
    plot(rays.RayYRList(rays.RayIndList(j):rays.RayEndList(j),1),...
         rays.RayYRList(rays.RayIndList(j):rays.RayEndList(j),2),...
         'b-','lineWidth',2)
    % plot ic of the rays
    plot(rays.RayYRList(rays.RayIndList(j),1),...
         rays.RayYRList(rays.RayIndList(j),2),'bo')
end


%% Plot k and group vel vectors along a ray
j=1;
rayx=rays.RayYRList(rays.RayIndList(j):rays.RayEndList(j),1);
rayz=rays.RayYRList(rays.RayIndList(j):rays.RayEndList(j),2);
raykx=rays.RayYRList(rays.RayIndList(j):rays.RayEndList(j),3);
raykz=rays.RayYRList(rays.RayIndList(j):rays.RayEndList(j),4);
time =rays.RayTRList(rays.RayIndList(j):rays.RayEndList(j));
dt = diff(time);
xdot = diff(rayx)./dt;
zdot = diff(rayz)./dt;

%%
figure(1)
plot(time(1:end-1),xdot,time(1:end-1),zdot)
 
 %%
figure(1)
plot(time(1:end-1),((xdot.*raykx(2:end)+zdot.*raykz(2:end))  ));

%%   
    
figure(2)
clf 
%plasma=plotFlux(plasma);
hold on
plot(rayx,rayz,'b-','lineWidth',2)
% plot ic of the rays
plot(rayx(1),rayz(1),'bo')
% quiver(rayx(1:5:end),rayz(1:5:end),...
%        raykx(1:5:end),raykz(1:5:end),'ShowArrowHead','off')
 quiver(rayx(1:5:end-1),rayz(1:5:end-1),...
        xdot(1:5:end),zdot(1:5:end),'ShowArrowHead','off')


function rays = initRays()     % Raytracing variables
    rays =struct(...     % 
    'NRAY',[], ...   %   number of rays to propagate
    'MON',[], ...    %   unsorted monitors accumulated during the evolution
    'MONNorm',[], ...%   normalizing factors
    'TYPE',[], ...   %   calculation (Trj=trajectory,Con=conversion,Amp=amplitude)
    'aThres',[], ... %   amplitude threshold below which rays are removed
    'convert',[], ...%   list of indices for rays to convert
    'hp',[], ...     %   5 steps history for time derivatives, value polarizations
    'ht',[], ...     %   5 steps history for time derivatives, value of t
    'hy',[], ...     %   5 steps history for time derivatives, value of y=(x,k,S,A,T)
    'inKspace',[],...%   flag for evolution of each ray (0=config space, 1=k-space)
    'monctc',[], ... %   last recorded caustics monitor
    'moncnv',[], ... %   last recorded conversion monitor
    'moncnv_der',[],... 
    'plt',[], ...    %   components of ODE to plot
    'stp',[], ...    %   index to monitors along the ray
    'caustic',[],... %   list of indices for rays to perform Maslov transform
    'odeDim',[], ... %   dimension of ODE vector
    'odeOptions',[],... %   controlling ODE solver
    'NE',[], ...        %   event ray number (which ray converts)
    'TR',[], ...     %   solution times (rays' trajectories)
    'YR',[], ...     %   solution components (rays' trajectories)
    'time',[], ...   %   current time
    'timespan',[], ... % duration of complete calculation
    'timeintv',[], ... % duration until next stop
    'tspan',[], ...  %   time interval tspan = [time time+timeintv]
    'y',[],   ...    %   ODE components
    'kray0',[], ...  %   ray pencil initial wave vectors in (r,phi,z)
    'sray0',[], ...  %   ray pencil initial normalized radii
    'thray0',[],...  %   ray pencil initial poloidal angles
    'yalf0',[],...
    'eta2est',[],...
    'adjcnv',[],...
    'adjctc',[],...
    'adjinit',1,...
    'RayICList',[],... % holds initial conditions for rays that need to be traced
    'RayIniTimeList',[],... % holds initial time for the rays.
    'RayTRList',[],... % holds time var from the raytracing for each ray
    'RayYRList',[],... % holds the trajectory for each ray
    'RayIndList', [],...% holds the index where the jth ray's data starts in RayTRList
    'RayPhaseList',[],...
    'RayInKList',[],...
    'RayConvList',[],... % holds data about mode converted rays
    'monErrAbort',0,...
    'initialstep',[],...
    'zstList',[]); % list of saddle points computed during mode conversion calculations

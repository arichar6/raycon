function plasma = initPlasma()
    plasma =struct(...    % A struct which will contain the system config variables
    'EQ',[], ...     %   equilibrium (either 'Solovev' or a filename)
    'MODEL',[], ...  %   dispersion equation model (msw1x1, cld2x2, cld3x3, etc)
    'PROBL',[], ...  %   problem (tok=tokamak, oth=Steve define your own)
    'acharge',[], ...%   atomic charge of species (e.g. -1 for electron)
    'amass',[], ...  %   atomic mass of species (e.g. 4 for helium)
    'b0',[], ...     %   magnetic field on axis
    'elong',[], ...  %   plasma elongation
    'freq',[], ...   %   antenna linear frequency
    'iaspr',[], ...  %   plasma inverse aspect ratio
    'psin',[], ...   %   magnetic flux at the plasma edge for normalization
    'q0',[], ...     %   plasma safety factor on axis
    'r0',[], ...     %   plasma major radius
    'n0',[], ...     %   densities on axis
    'na',[], ...     %   densities profile factor in n=n0*(1-na*s^2)^nb
    'nb',[], ...     %   densities profile factor in n=n0*(1-na*s^2)^nb
    't0',[], ...     %   temperatures on axis
    'ta',[], ...     %   temperatures profile factor in t=t0*(1-ta*s^2)^tb
    'tb',[], ...     %   temperatures profile factor in t=t0*(1-ta*s^2)^tb
    'kant',[], ...   %   antenna wave vector in (r,phi,z)
    'sant',[], ...   %   antenna normalized radius
    'thant',[], ...  %   antenna lower and upper poloidal angles
    'omega',[], ...  %   antenna circular frequency
    'r',[],...
    'z',[],...
    'Rho',[],...
    'Theta',[],...
    's',[], ...      %   normalized radial mesh
    'theta',[], ...  %   poloidal mesh
    'depo',[]);   %   radial deposition profiles per species


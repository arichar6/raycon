function cnst = initCnst()
    cnst = struct(...    %
     'c',[], ...      %   speed of light in vacuum [m/s]
     'e',[], ...      %   electron charge [C]
     'eps0',[], ...   %   permittivity of vacuum [F/m]
     'mp',[]) ;        %   proton mass [kg]
 
    cnst.c        = 2.9979E+8;                      %   speed of light in vacuum [m/s]
    cnst.e        = 1.6022E-19;                     %   electron charge [C]
    cnst.mp       = 1.6726E-27;                     %   proton mass [kg]
    cnst.eps0     = 8.8542E-12;                     %   permittivity of vacuum [F/m]

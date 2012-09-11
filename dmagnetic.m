function [dbdr,dbdz,dbdr2,dbdrz,dbdz2] = dmagnetic(rho,theta)

% NB, the input is required to be scalar.

global plasma
%
% ----- Coordinates
r= rho*cos(theta)+plasma.r0;
z=-rho*sin(theta);

%
% ----- Size of the step used to estimate the derivatives
dr = 1E-8;
dz = 1E-8;

[R,Z] = meshgrid((r-dr):dr:(r+dr), (z-dz):dz:(z+dz));
Rho   = sqrt((R-plasma.r0).^2+Z.^2);
Theta = atan2(plasma.r0-R,-Z)+pi/2;

%
% ----- Radial magnetic flux variable
if strcmp(plasma.EQ,'Solovev')
  [sflx,dsdr,dsdz] = solovev(Rho,Theta,plasma.r0,plasma.iaspr,plasma.elong,0);
  dpds =2.*sflx*plasma.psin;  
  t0=plasma.b0.*plasma.r0;
end

%
% ----- Magnetic field amplitude
h11   = dsdr.^2+dsdz.^2;
gp2   = dpds.^2.*h11;
bpol2 = gp2./(R.^2);               
btor2 = (t0./R).^2;
b     = sqrt(btor2 + bpol2);

%
% ----- Numerical estimate of the derivatives
dbdr  = (b(2,3)-b(2,1))/(2*dr);
dbdz  = (b(3,2)-b(1,2))/(2*dz);
dbdr2 = (b(2,3)-2*b(2,2)+b(2,1))/(dr^2);
dbdz2 = (b(3,2)-2*b(2,2)+b(1,2))/(dz^2);
dbdrz = (b(3,3)-b(3,1)-b(1,3)+b(1,1))/(4*dr*dz);

return
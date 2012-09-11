function [rho,r,z] = mapFlux(s,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MAPFLUX -- Map from flux coordinates to real space
%             Part of the RAYcON package
%
%  A. JAUN, Alfven Laboratory, KTH, 100 44 Stockholm, Sweden
%  A.N. KAUFMAN, Lawrence Berkeley Laboratory, Berkeley, CA 94720, USA
%  E.R. TRACY, College of William & Mary, Williamsburg, VA 23187-8795, USA
%
%  (C) Version 7.0,  14-Aug-2006. All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global plasma rays sys cnst
%
if strcmp(plasma.EQ,'Solovev')                     % Analytical equilibrium
  rhoa=plasma.elong*plasma.iaspr*plasma.r0;        % Estimate radius
  rho =zeros(size(theta)); %der=1;
  for j=1:size(theta,2);
    rho(j)=fzero('solovev',[0 1.5*rhoa], ...
                 optimset('TolX',1e-4,'Display','off'), ...
                 theta(j),plasma.r0,plasma.iaspr,plasma.elong,s(j));
  end;
  r= rho.*cos(theta)+plasma.r0;
  z=-rho.*sin(theta);
end;




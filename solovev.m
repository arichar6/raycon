function [varargout] = solovev(rho,theta,r0,eps,elong,sflxa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SOLOVEV -- Solovev equilibrium in normalized flux s=sqrt(psi/psia)
%             Part of the RAYcON package
%
%  A. JAUN, Alfven Laboratory, KTH, 100 44 Stockholm, Sweden
%  A.N. KAUFMAN, Lawrence Berkeley Laboratory, Berkeley, CA 94720, USA
%  E.R. TRACY, College of William & Mary, Williamsburg, VA 23187-8795, USA
%
%  (C) Version 7.0,  14-Aug-2006. All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global plasma

%
% ----- Cylindrical coordinates and aliases
%r= rho.*cos(theta)+r0; rsq=r.^2;  psin=1.;   Esq=elong^2;  r0sq=r0^2;
r= rho.*cos(theta)+r0; rsq=r.^2; psin=plasma.psin; Esq=elong^2;  r0sq=r0^2;
z=-rho.*sin(theta);    zsq=z.^2;  fac =psin./(eps*r0sq)^2;
%
% ----- Flux function
psi=fac*(rsq.*zsq/Esq+0.25*(rsq-r0sq).^2); 
sflx=sqrt(psi/psin)-sflxa;   % NB: sflxa is used to find the point where sflx=sflxa
sflx2=sflx.^2;
if (nargout==1) 
  [varargout{1:1}]=deal(sflx);
  return
end

% ----- 1st derivatives
dpdr  =fac*r.*(2*zsq./Esq + (rsq-r0sq));     dpdrsq=dpdr.^2;
dpdz  =fac*z.*(2*rsq./Esq);                  dpdzsq=dpdz.^2;
dpds = 2.*sflx*psin; dsdp =1./dpds;
dsdr = dsdp.*dpdr;
dsdz = dsdp.*dpdz;

if (nargout==3) 
  [varargout{1:3}]=deal(sflx,dsdr,dsdz);
  return
end

% ----- 2nd derivatives
dpdr2 =fac*(2*zsq./Esq + 3*rsq - r0sq);
dpdrz =fac*(4*r.*z./Esq);
dpdz2 =fac*(2*rsq./Esq);
dsdr2= dsdp.*(dpdr2-dsdp./sflx.*dpdrsq);
dsdrz= dsdp.*(dpdrz-dsdp./sflx.*dpdr.*dpdz);
dsdz2= dsdp.*(dpdz2-dsdp./sflx.*dpdzsq);

if (nargout==6) 
  [varargout{1:6}]=deal(sflx,dsdr,dsdz,dsdr2,dsdrz,dsdz2);
  return
end

% ----- 3rd derivatives
dpdr3 =fac*6*r;
dpdr2z=fac*4*z./Esq;
dpdrz2=fac*4*r./Esq;
dpdz3 =zeros(size(psi));

dLNsdp= dsdp./sflx;
dsdrp =-dLNsdp.*dsdp.*dpdr; dLNsdrp=dLNsdp.*(dsdrp./dsdp-dsdr./sflx);
dsdzp =-dLNsdp.*dsdp.*dpdz; dLNsdzp=dLNsdp.*(dsdzp./dsdp-dsdz./sflx);

dsdr3 =dsdr2.*(dsdrp./dsdp)+dsdp.*(dpdr3 -dLNsdrp.*dpdrsq-2*dLNsdp.*dpdr.*dpdr2);
dsdr2z=dsdr2.*(dsdzp./dsdp)+dsdp.*(dpdr2z-dLNsdzp.*dpdrsq-2*dLNsdp.*dpdr.*dpdrz);
dsdrz2=dsdz2.*(dsdrp./dsdp)+dsdp.*(dpdrz2-dLNsdrp.*dpdzsq-2*dLNsdp.*dpdz.*dpdrz);
dsdz3 =dsdz2.*(dsdzp./dsdp)+dsdp.*(dpdz3 -dLNsdzp.*dpdzsq-2*dLNsdp.*dpdz.*dpdz2);

[varargout{1:10}]= ...
   deal(sflx,dsdr,dsdz,dsdr2,dsdrz,dsdz2,dsdr3,dsdr2z,dsdrz2,dsdz3);
return;

function [varargout] = magnetic(rho,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MAGNETIC --  Magnetic field and topology.
%               Part of the RAYcON package
%
%  A. JAUN, Alfven Laboratory, KTH, 100 44 Stockholm, Sweden
%  A.N. KAUFMAN, Lawrence Berkeley Laboratory, Berkeley, CA 94720, USA
%  E.R. TRACY, College of William & Mary, Williamsburg, VA 23187-8795, USA
%
%  (C) Version 7.0,  14-Aug-2006. All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global plasma
%
% ----- Coordinates
cost=cos(theta); sint=sin(theta); rho2=rho.^2;
r= rho.*cost+plasma.r0; r2=r.^2;
z=-rho.*sint;    z2=z.^2;
%
% ----- Radial magnetic flux variable
if strcmp(plasma.EQ,'Solovev')
  [sflx,dsdr,dsdz,dsdr2,dsdrz,dsdz2,dsdr3,dsdr2z,dsdrz2,dsdz3]= ...
     solovev(rho,theta,plasma.r0,plasma.iaspr,plasma.elong,0);
  psi=plasma.psin.*sflx.^2; 
  dpds =2.*sflx*plasma.psin;  
  dsdp = 1./dpds;
  dpds2=2.*plasma.psin;
  dsdp2=-2*dsdp/sflx;
  t0=plasma.b0.*plasma.r0; dt0ds=zeros(size(plasma.r0));
  dpds2=2*plasma.psin;

end;
%
% ----- Poloidal angle
dtdr =-sint./rho;
dtdz =-cost./rho;     
dtdr2= 2*sint.*cost./rho2;
dtdrz= (cost.^2-sint.^2)./rho2;
dtdz2=-2*sint.*cost./rho2;
%
% ----- Metric coefficients (1st order)
jac   = r./(dsdz.*dtdr-dsdr.*dtdz);
drds  =-jac./r.*dtdz;  dzds = jac./r.*dtdr; 
drdt  = jac./r.*dsdz;  dzdt =-jac./r.*dsdr;
h11   = dsdr.^2+dsdz.^2;
dh11dr= 2*(dsdr.*dsdr2 + dsdz.*dsdrz); 
dh11dz= 2*(dsdr.*dsdrz + dsdz.*dsdz2);
dh11ds= dh11dr.*drds + dh11dz.*dzds;
dh11dt= dh11dr.*drdt + dh11dz.*dzdt;
gp2   = dpds.^2.*h11;
dgp2ds= dpds.^2.*(dh11ds+2*h11./sflx);
dgp2dt= dpds.^2.* dh11dt;
%
% ----- Magnetic field amplitude (1st derivatives)
a2    = t0.^2+gp2;
da2ds = 2*t0.*dt0ds + dgp2ds;
da2dt =               dgp2dt;
bpol2 = gp2./r2;               
btor2=(t0./r).^2;
b     = sqrt(btor2 + bpol2);
db1ds = da2ds./(2*b.*r2); db2ds=-b.*drds./r; dbds=db1ds+db2ds;
db1dt = da2dt./(2*b.*r2); db2dt=-b.*drdt./r; dbdt=db1dt+db2dt;
bp=sqrt(bpol2./btor2);
%
% ----- Magnetic field orientation (1st derivatives)
gs  = sqrt(h11);
ener= dsdr./gs;
enef= zeros(size(b));
enez= dsdz./gs;
eber= t0.*dsdz./(r.*b.*gs);
ebef=-dpds.*gs./(r.*b);
ebez=-t0.*dsdr./(r.*b.*gs);
eper= dpds.*dsdz./(r.*b);
epef= t0./(r.*b);
epez=-dpds.*dsdr./(r.*b);
%
if (nargout==18) 
 [varargout{1:18}]=deal(b,dbds,dbdt,bp,sflx,dsdr,dsdz,dtdr,dtdz, ...
                        ener,enef,enez,eber,ebef,ebez,eper,epef,epez);
 return;
end;

% dsdrs = dsdr2.*drds+dsdrz.*dzds;
dsdrt = dsdr2.*drdt+dsdrz.*dzdt;
% dsdzs = dsdrz.*drds+dsdz2.*dzds;
dsdzt = dsdrz.*drdt+dsdz2.*dzdt;
dtdrs = dtdr2.*drds+dtdrz.*dzds;
dtdrt = dtdr2.*drdt+dtdrz.*dzdt;
dtdzs = dtdrz.*drds+dtdz2.*dzds;
dtdzt = dtdrz.*drdt+dtdz2.*dzdt;

% ----- Metric coefficients (2nd order) %Olga opposite sign after jac.*
djacdr=jac./r.*(1-jac.*(dsdrz.*dtdr+dsdz.*dtdr2-dsdr2.*dtdz-dsdr.*dtdrz)); 
djacdz=jac./r.*( -jac.*(dsdz2.*dtdr+dsdz.*dtdrz-dsdrz.*dtdz-dsdr.*dtdz2));
djacds=djacdr.*drds+djacdz.*dzds;    djacdt=djacdr.*drdt+djacdz.*dzdt;

drds2 =-(djacds.*dtdz-jac.*drds.*dtdz./r+jac.*dtdzs)./r;
dzds2 = (djacds.*dtdr-jac.*drds.*dtdr./r+jac.*dtdrs)./r;
drdst =-(djacdt.*dtdz-jac.*drdt.*dtdz./r+jac.*dtdzt)./r;
dzdst = (djacdt.*dtdr-jac.*drdt.*dtdr./r+jac.*dtdrt)./r;
drdt2 = (djacdt.*dsdz-jac.*drdt.*dsdz./r+jac.*dsdzt)./r;
dzdt2 =-(djacdt.*dsdr-jac.*drdt.*dsdr./r+jac.*dsdrt)./r;

dh11dr2=2*(dsdr2.*dsdr2+dsdr.*dsdr3 +dsdrz.*dsdrz+dsdz.*dsdr2z);
dh11drz=2*(dsdrz.*dsdr2+dsdr.*dsdr2z+dsdz2.*dsdrz+dsdz.*dsdrz2);
dh11dz2=2*(dsdrz.*dsdrz+dsdr.*dsdrz2+dsdz2.*dsdz2+dsdz.*dsdz3 );
% dh11drs=dh11dr2.*drds+dh11drz.*dzds;  dh11drt=dh11dr2.*drdt+dh11drz.*dzdt;
% dh11dzs=dh11drz.*drds+dh11dz2.*dzds;  dh11dzt=dh11drz.*drdt+dh11dz2.*dzdt;
% dh11ds2=dh11drs.*drds+dh11dr.*drds2 +dh11dzs.*dzds+dh11dz.*dzds2;
% dh11dst=dh11drs.*drdt+dh11dr.*drdst +dh11dzs.*dzdt+dh11dz.*dzdst;
% dh11dt2=dh11drt.*drdt+dh11dr.*drdt2 +dh11dzt.*dzdt+dh11dz.*dzdt2;

% dgp2ds2=dpds.^2.*(dh11ds2+2*dh11ds./sflx+2*h11./sflx.^2)+2*dpds2*dgp2ds./dpds;
% dgp2dst=dpds.^2.*(dh11dst+2*dh11dt./sflx);
% dgp2dt2=dpds.^2.* dh11dt2 ;

% ----- Magnetic field amplitude (2nd derivatives)
% da2ds2=dgp2ds2;
% da2dst=dgp2dst;
% da2dt2=dgp2dt2;

% db1ds2= db1ds.*(da2ds2./da2ds -(dbds.*r2 +2*b.*r.*drds)./(b.*r2));
% db1dst= db1ds.*(da2dst./da2ds -(dbdt.*r2 +2*b.*r.*drdt)./(b.*r2));
% db1dt2= db1dt.*(da2dt2./da2dt -(dbdt.*r2 +2*b.*r.*drdt)./(b.*r2));
% db2ds2= db2ds.*(dbds./b +drds2./drds -drds./r); 
% db2dst= db2ds.*(dbdt./b +drdst./drds -drdt./r); 
% db2dt2= db2dt.*(dbdt./b +drdt2./drdt -drdt./r); 
% dbds2=db1ds2+db2ds2;
% dbdst=db1dst+db2dst;
% dbdt2=db1dt2+db2dt2;

% dbdr2 = dsdr2.*dbds + dbds2.*(dsdr).^2 + 2*dtdr.*dsdr.*dbdst ...
%         + dbdt2.*(dtdr).^2 + dtdr2.*dbdt;
% dbdz2 = dsdz2.*dbds + dbds2.*(dsdz).^2 + 2*dtdz.*dsdz.*dbdst ...
%         + dbdt2.*(dtdz).^2 + dtdz2.*dbdt;
% dbdrz = dsdrz.*dbds + dtdrz.*dbdt + dsdr.*(dsdz.*dbds2 + dtdz.*dbdst) ...
%         + dtdr.*(dsdz.*dbdst + dtdz.*dbdt2);

if numel(rho)==1
    [dbdr,dbdz,dbdr2,dbdrz,dbdz2] = dmagnetic(rho,theta);
else % numel(rho) > 1
    dbdr = zeros(size(rho,1),size(rho,2));
    dbdz = zeros(size(rho,1),size(rho,2));
    dbdr2 = zeros(size(rho,1),size(rho,2));
    dbdrz = zeros(size(rho,1),size(rho,2));
    dbdz2 = zeros(size(rho,1),size(rho,2));
    for j=1:size(rho,1)
        for k=1:size(rho,2)
        [ndbdr,ndbdz,ndbdr2,ndbdrz,ndbdz2] = dmagnetic(rho(j,k),theta(j,k));
        dbdr(j,k) = ndbdr;
        dbdz(j,k) = ndbdz;
        dbdr2(j,k) = ndbdr2;
        dbdrz(j,k) = ndbdrz;
        dbdz2(j,k) = ndbdz2;
        end
    end
end

dbds2 = drds2.*dbdr + dzds2.*dbdz + dbdr2.*drds.^2 + dbdz2.*dzds.^2 + 2*drds.*dzds.*dbdrz;
dbdt2 = drdt2.*dbdr + dzdt2.*dbdz + dbdr2.*drdt.^2 + dbdz2.*dzdt.^2 + 2*drdt.*dzdt.*dbdrz;
dbdst = drdst.*dbdr + dzdst.*dbdz + drds.*drdt.*dbdr2 + drds.*dzdt.*dbdrz ...
        + dzds.*drdt.*dbdrz + dzds.*dzdt.*dbdz2;

if (nargout==31) 
 [varargout{1:31}]= deal(b,dbds,dbdt,bp,sflx,dsdr,dsdz,dtdr,dtdz, ...
                         ener,enef,enez,eber,ebef,ebez,eper,epef,epez, ...
                         dbds2,dbdst,dbdt2, ...
                         dsdr2,dsdrz,dsdz2,dsdr3,dsdr2z,dsdrz2,dsdz3, ...
                         dtdr2,dtdrz,dtdz2);
 return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Added by Steve R.  Field curvature terms.           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% zeta = (r*b)^(-1);
zeta = (r.*b)^(-1);
dzetadr = -zeta.*(1./r +dbdr./b);
dzetadz = -zeta.*dbdz./b;

% mag = 1/sqrt(h11) = 1/gs
mag = 1./gs;
dmagdr = -(1/2)*(mag.^3).*dh11dr;
dmagdz = -(1/2)*(mag.^3).*dh11dz;

% Derivatives
denerdr = (1./gs).*dsdr2 + dsdr.*dmagdr;
denefdr = zeros(size(enef));
denezdr = (1./gs).*dsdrz + dsdz.*dmagdr;

denerdz = (1./gs).*dsdrz + dsdr.*dmagdz;
denefdz = zeros(size(enef));
denezdz = (1./gs).*dsdz2 + dsdz.*dmagdz;

deberdr = t0*(dsdrz.*zeta./gs + dsdz.*dzetadr./gs + dsdz.*zeta.*dmagdr);
debefdr = -(dsdr.*dpds2.*zeta.*gs + dpds.*dzetadr.*gs +...
            dpds.*zeta.*(dsdr.*dsdr2+dsdz.*dsdrz)./gs);
debezdr = -t0*(dsdr2.*zeta./gs + dsdr.*dzetadr./gs + dsdr.*zeta.*dmagdr);

deberdz = t0*(dsdz2.*zeta./gs + dsdz.*dzetadz./gs + dsdz.*zeta.*dmagdz);
debefdz = -(dsdz.*dpds2.*zeta.*gs + dpds.*dzetadz.*gs +...
            dpds.*zeta.*(dsdz.*dsdz2+dsdr.*dsdrz)./gs);
debezdz = -t0*(dsdrz.*zeta./gs + dsdr.*dzetadz./gs + dsdr.*zeta.*dmagdz);

deperdr = dsdr.*dpds2.*dsdz.*zeta + dpds.*dsdrz.*zeta + dpds.*dsdz.*dzetadr;
depefdr = t0*dzetadr;
depezdr =-dsdr.*dpds2.*dsdr.*zeta - dpds.*dsdr2.*zeta - dpds.*dsdr.*dzetadr;

deperdz = dsdz.*dpds2.*dsdz.*zeta + dpds.*dsdz2.*zeta + dpds.*dsdz.*dzetadz;
depefdz = t0*dzetadz;
depezdz =-dsdz.*dpds2.*dsdr.*zeta - dpds.*dsdrz.*zeta - dpds.*dsdr.*dzetadz;

if (nargout==49) 
    [varargout{1:49}]= deal(b,dbds,dbdt,bp,sflx,dsdr,dsdz,dtdr,dtdz, ...
                         ener,enef,enez,eber,ebef,ebez,eper,epef,epez, ...
                         dbds2,dbdst,dbdt2, ...
                         dsdr2,dsdrz,dsdz2,dsdr3,dsdr2z,dsdrz2,dsdz3, ...
                         dtdr2,dtdrz,dtdz2,...
                         denerdr, denefdr, denezdr,...
                         denerdz, denefdz, denezdz,...
                         deberdr, debefdr, debezdr,...
                         deberdz, debefdz, debezdz,...
                         deperdr, depefdr, depezdr,...
                         deperdz, depefdz, depezdz...
                         );
    return
end

% Need _2nd_ derivatives of the ener... elements!

% First calculate the derivatives of some of the terms in the expressions

dzetadr2 = -dzetadr.*(1./r + dbdr./b) ...
           -zeta.*(-1./r.^2 -(dbdr./b).^2 + dbdr2./b);
dzetadrz = -dzetadr.*dbdz./b -zeta.*(dbdrz./b - dbdr.*dbdz./(b.^2));
dzetadz2 = -dzetadz.*dbdz./b -zeta.*(dbdz2./b - (dbdz./b).^2);

% dzetadr2 = -(1./r).*(dbdr2./b.^2 -2*dbdr.^2./b.^3) + 2*dbdr./(r.^2.*b.^2)...
%             +2/(b.*r.^3);
% dzetadrz = -(1./r).*(dbdrz./b.^2 -2*dbdr.*dbdz./b.^3) + dbdz./(r.^2.*b.^2);
% dzetadz2 = -(1./r).*(dbdz2./b.^2 -2*dbdz.^2./b.^3);

dmagdr2 = (3/4)*(mag.^5).*(dh11dr.^2) - (1/2)*dh11dr2.*mag.^3;
dmagdrz = (3/4)*(mag.^5).*(dh11dr.*dh11dz) - (1/2)*dh11drz.*mag.^3;
dmagdz2 = (3/4)*(mag.^5).*(dh11dz.^2) - (1/2)*dh11dz2.*mag.^3;

dgsdr = (1/2)*mag.*dh11dr;
dgsdz = (1/2)*mag.*dh11dz;
dgsdr2 = (1/2)*(dmagdr.*dh11dr + mag.*dh11dr2);
dgsdrz = (1/2)*(dmagdr.*dh11dz + mag.*dh11drz);
dgsdz2 = (1/2)*(dmagdz.*dh11dz + mag.*dh11dz2);

dpdsr = dsdr.*dpds2;
dpdsz = dsdz.*dpds2;
dpdsr2 = dsdr2.*dpds2;
dpdsrz = dsdrz.*dpds2;
dpdsz2 = dsdz2.*dpds2;

% ener
denerdr2 = mag.*dsdr3 + 2*dmagdr.*dsdr2 + dsdr.*dmagdr2;
denerdrz = mag.*dsdr2z + dmagdz.*dsdr2 + dmagdr.*dsdrz + dsdr.*dmagdrz;
denerdz2 = mag.*dsdrz2 + 2*dmagdz.*dsdrz + dsdr.*dmagdz2;

% enef
denefdr2 = zeros(size(enef));
denefdrz = zeros(size(enef));
denefdz2 = zeros(size(enef));

% enez
denezdr2 = dmagdr2.*dsdz + 2*dmagdr.*dsdrz + mag.*dsdr2z;
denezdrz = dmagdrz.*dsdz + dmagdr.*dsdz2 + dmagdz.*dsdrz + mag.*dsdrz2;
denezdz2 = dmagdz2.*dsdz + 2*dmagdz.*dsdz2 + mag.*dsdz3;

% eber
% deberdij = t0*(dsdzij.*zeta.*mag + dsdzj.*dzetadi.*mag + dsdzj.*zeta.*dmagdi ...
%              + dsdzi.*dzetadj.*mag + dsdz.*dzetadij.*mag + dsdz.*dzetadj.*dmagdi ...
%              + dsdzi.*zeta.*dmagdj + dsdz.*dzetadi.*dmagdj + dsdz.*zeta.*dmagdij);
deberdr2 = t0*(dsdr2z.*zeta.*mag + dsdrz.*dzetadr.*mag + dsdrz.*zeta.*dmagdr ...
             + dsdrz.*dzetadr.*mag + dsdz.*dzetadr2.*mag + dsdz.*dzetadr.*dmagdr ...
             + dsdrz.*zeta.*dmagdr + dsdz.*dzetadr.*dmagdr + dsdz.*zeta.*dmagdr2);
deberdrz = t0*(dsdrz2.*zeta.*mag + dsdz2.*dzetadr.*mag + dsdz2.*zeta.*dmagdr ...
             + dsdrz.*dzetadz.*mag + dsdz.*dzetadrz.*mag + dsdz.*dzetadz.*dmagdr ...
             + dsdrz.*zeta.*dmagdz + dsdz.*dzetadr.*dmagdz + dsdz.*zeta.*dmagdrz);
deberdz2 = t0*(dsdz3.*zeta.*mag + dsdz2.*dzetadz.*mag + dsdz2.*zeta.*dmagdz ...
             + dsdz2.*dzetadz.*mag + dsdz.*dzetadz2.*mag + dsdz.*dzetadz.*dmagdz ...
             + dsdz2.*zeta.*dmagdz + dsdz.*dzetadz.*dmagdz + dsdz.*zeta.*dmagdz2);

% ebef
% debefdij = -(dpdsij.*zeta.*gs + dpdsj.*dzetadi.*gs + dpdsj.*zeta.*dgsdi ...
%              + dpdsi.*dzetadj.*gs + dpds.*dzetadij.*gs + dpds.*dzetadj.*dgsdi ...
%              + dpdsi.*zeta.*dgsdj + dpds.*dzetadi.*dgsdj + dpds.*zeta.*dgsdij);
debefdr2 = -(dpdsr2.*zeta.*gs + dpdsr.*dzetadr.*gs + dpdsr.*zeta.*dgsdr ...
             + dpdsr.*dzetadr.*gs + dpds.*dzetadr2.*gs + dpds.*dzetadr.*dgsdr ...
             + dpdsr.*zeta.*dgsdr + dpds.*dzetadr.*dgsdr + dpds.*zeta.*dgsdr2);
debefdrz = -(dpdsrz.*zeta.*gs + dpdsz.*dzetadr.*gs + dpdsz.*zeta.*dgsdr ...
             + dpdsr.*dzetadz.*gs + dpds.*dzetadrz.*gs + dpds.*dzetadz.*dgsdr ...
             + dpdsr.*zeta.*dgsdz + dpds.*dzetadr.*dgsdz + dpds.*zeta.*dgsdrz);
debefdz2 = -(dpdsz2.*zeta.*gs + dpdsz.*dzetadz.*gs + dpdsz.*zeta.*dgsdz ...
             + dpdsz.*dzetadz.*gs + dpds.*dzetadz2.*gs + dpds.*dzetadz.*dgsdz ...
             + dpdsz.*zeta.*dgsdz + dpds.*dzetadz.*dgsdz + dpds.*zeta.*dgsdz2);

% ebez
% debezdij = -t0*(dsdrij.*zeta.*mag + dsdrj.*dzetadi.*mag + dsdrj.*zeta.*dmagdi ...
%              + dsdri.*dzetadj.*mag + dsdr.*dzetadij.*mag + dsdr.*dzetadj.*dmagdi ...
%              + dsdri.*zeta.*dmagdj + dsdr.*dzetadi.*dmagdj + dsdr.*zeta.*dmagdij);
debezdr2 = -t0*(dsdr3.*zeta.*mag + dsdr2.*dzetadr.*mag + dsdr2.*zeta.*dmagdr ...
             + dsdr2.*dzetadr.*mag + dsdr.*dzetadr2.*mag + dsdr.*dzetadr.*dmagdr ...
             + dsdr2.*zeta.*dmagdr + dsdr.*dzetadr.*dmagdr + dsdr.*zeta.*dmagdr2);
debezdrz = -t0*(dsdr2z.*zeta.*mag + dsdrz.*dzetadr.*mag + dsdrz.*zeta.*dmagdr ...
             + dsdr2.*dzetadz.*mag + dsdr.*dzetadrz.*mag + dsdr.*dzetadz.*dmagdr ...
             + dsdr2.*zeta.*dmagdz + dsdr.*dzetadr.*dmagdz + dsdr.*zeta.*dmagdrz);
debezdz2 = -t0*(dsdrz2.*zeta.*mag + dsdrz.*dzetadz.*mag + dsdrz.*zeta.*dmagdz ...
             + dsdrz.*dzetadz.*mag + dsdr.*dzetadz2.*mag + dsdr.*dzetadz.*dmagdz ...
             + dsdrz.*zeta.*dmagdz + dsdr.*dzetadz.*dmagdz + dsdr.*zeta.*dmagdz2);

% eper
% deperdij = dsdzij.*zeta.*dpds + dsdzj.*dzetadi.*dpds + dsdzj.*zeta.*dpdsi ...
%          + dsdzi.*dzetadj.*dpds + dsdz.*dzetadij.*dpds + dsdz.*dzetadj.*dpdsi ...
%          + dsdzi.*zeta.*dpdsj + dsdz.*dzetadi.*dpdsj + dsdz.*zeta.*dpdsij;
deperdr2 = dsdr2z.*zeta.*dpds + dsdrz.*dzetadr.*dpds + dsdrz.*zeta.*dpdsr ...
         + dsdrz.*dzetadr.*dpds + dsdz.*dzetadr2.*dpds + dsdz.*dzetadr.*dpdsr ...
         + dsdrz.*zeta.*dpdsr + dsdz.*dzetadr.*dpdsr + dsdz.*zeta.*dpdsr2;
deperdrz = dsdrz2.*zeta.*dpds + dsdz2.*dzetadr.*dpds + dsdz2.*zeta.*dpdsr ...
         + dsdrz.*dzetadz.*dpds + dsdz.*dzetadrz.*dpds + dsdz.*dzetadz.*dpdsr ...
         + dsdrz.*zeta.*dpdsz + dsdz.*dzetadr.*dpdsz + dsdz.*zeta.*dpdsrz;
deperdz2 = dsdz3.*zeta.*dpds + dsdz2.*dzetadz.*dpds + dsdz2.*zeta.*dpdsz ...
         + dsdz2.*dzetadz.*dpds + dsdz.*dzetadz2.*dpds + dsdz.*dzetadz.*dpdsz ...
         + dsdz2.*zeta.*dpdsz + dsdz.*dzetadz.*dpdsz + dsdz.*zeta.*dpdsz2;

% epef
depefdr2 = t0*dzetadr2;
depefdrz = t0*dzetadrz;
depefdz2 = t0*dzetadz2;

% epez
depezdr2 = -(dsdr3.*zeta.*dpds + dsdr2.*dzetadr.*dpds + dsdr2.*zeta.*dpdsr ...
         + dsdr2.*dzetadr.*dpds + dsdr.*dzetadr2.*dpds + dsdr.*dzetadr.*dpdsr ...
         + dsdr2.*zeta.*dpdsr + dsdr.*dzetadr.*dpdsr + dsdr.*zeta.*dpdsr2);
depezdrz = -(dsdr2z.*zeta.*dpds + dsdrz.*dzetadr.*dpds + dsdrz.*zeta.*dpdsr ...
         + dsdr2.*dzetadz.*dpds + dsdr.*dzetadrz.*dpds + dsdr.*dzetadz.*dpdsr ...
         + dsdr2.*zeta.*dpdsz + dsdr.*dzetadr.*dpdsz + dsdr.*zeta.*dpdsrz);
depezdz2 = -(dsdrz2.*zeta.*dpds + dsdrz.*dzetadz.*dpds + dsdrz.*zeta.*dpdsz ...
         + dsdrz.*dzetadz.*dpds + dsdr.*dzetadz2.*dpds + dsdr.*dzetadz.*dpdsz ...
         + dsdrz.*zeta.*dpdsz + dsdr.*dzetadz.*dpdsz + dsdr.*zeta.*dpdsz2);


[varargout{1:76}]= deal(b,dbds,dbdt,bp,sflx,dsdr,dsdz,dtdr,dtdz, ...
                     ener,enef,enez,eber,ebef,ebez,eper,epef,epez, ...
                     dbds2,dbdst,dbdt2, ...
                     dsdr2,dsdrz,dsdz2,dsdr3,dsdr2z,dsdrz2,dsdz3, ...
                     dtdr2,dtdrz,dtdz2,...
                     denerdr, denefdr, denezdr,...
                     denerdz, denefdz, denezdz,...
                     deberdr, debefdr, debezdr,...
                     deberdz, debefdz, debezdz,...
                     deperdr, depefdr, depezdr,...
                     deperdz, depefdz, depezdz,... % end of first derivs
                     denerdr2, denerdrz, denerdz2,...
                     denefdr2, denefdrz, denefdz2,...
                     denezdr2, denezdrz, denezdz2,...
                     deberdr2, deberdrz, deberdz2,...
                     debefdr2, debefdrz, debefdz2,...
                     debezdr2, debezdrz, debezdz2,...
                     deperdr2, deperdrz, deperdz2,...
                     depefdr2, depefdrz, depefdz2,...
                     depezdr2, depezdrz, depezdz2...
                     );



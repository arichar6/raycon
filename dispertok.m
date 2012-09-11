function out= dispertok(var,yv,ydotv,yddotv,poldot,oper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DISPERSION --  Plasma dispersion characteristics
%                 Part of the RAYcON package
%
%  var     variable used to adjust e.g. yv(1) for non-linear root finding
%  yv      rhs of ODE
%  ydotv   1st derivative of yv
%  yddotv  2nd derivative of yv
%  poldot  1st derivative of polarization
%  oper    'Frq' cyclotron frequencies (+ion-hybrid if 3 species)
%          'Msw' Dispersion relation approximated 1x1
%          'Dsp' Disperison relation according to MODEL
%          'Pol' Polarization vector according to MODEL
%          'Ant' Antenna IC
%          'Trj' RHS for trajectory
%          'Amp' RHS for amplitute
%          'Mon' Monitoring along the trajectory
%          'Mch' Monitoring matching parameters
%          'Sdl' Hyperbola saddle point
%          'Cnv' Converted ray parameters
%          'Trs' Transmitted ray parameters
%          'Sgn' Sign of dt/dsigma, set by dUdom
%
%  ToDo    Optimize computing
%           - change Dij(c) vectors into D(c,i,j) matrix
%           - remove terms that are zero using emacs search/substitute
%
%  A. JAUN, Alfven Laboratory, KTH, 100 44 Stockholm, Sweden
%  A.N. KAUFMAN, Lawrence Berkeley Laboratory, Berkeley, CA 94720, USA
%  E.R. TRACY, College of William & Mary, Williamsburg, VA 23187-8795, USA
%
%  (C) Version 7.0,  14-Aug-2006. All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global cnst sys
global plasma rays

isCmplx=0;
%
% ----- Vectorized input quantities ------------------------------------------
%
dimytot=numel(yv);
dimy   =dimytot/rays.odeDim;
y    =reshape(yv    ,rays.odeDim,dimy)';
ydot =reshape(ydotv ,rays.odeDim,dimy)';
yddot=reshape(yddotv,rays.odeDim,dimy)';
%
% ----- Correction 'var' is col vector to find zeros of dispersion function
yin=y+var*ydot;
%
% ----- Adjusted position & wave vector in cylindrical coordinates -----------
zro=zeros(size(yin(:,1))); one=ones(size(zro));
if (rays.odeDim==4||rays.odeDim>=9)
 rr=yin(:,1); zz=yin(:,2); kr=yin(:,3); kz=yin(:,4);
 kf=one.*plasma.kant(2);                            % fi=omega*tim/(kf*rr);
end
if(rays.odeDim>=9)
 dWr2=yin(:,5); dWrz=yin(:,6); dWz2=yin(:,7); lnE2 =yin(:,8); ph=yin(:,9);
end
rho=sqrt((rr-plasma.r0).^2+zz.^2);
theta=atan2(plasma.r0-rr,-zz)+pi/2;
%
% ----- Delimit calculation of higher order derivatives ----------------------
eval1st = ~(strcmp(oper,'Frq')|strcmp(oper,'Msw')|...
            strcmp(oper,'Dsp')|strcmp(oper,'Pol'));
eval2nd = ~(strcmp(oper,'Trj')|strcmp(oper,'Ant'))  & eval1st;
evalAll = strcmp(oper,'Amp');
%
%
% ===== Local plasma parameters ==============================================
%
% ----- Magnetic field and topology ------------------------------------------
if (~eval2nd)
  [b,dbds,dbdt,bp,sflx,dsdr,dsdz,dtdr,dtdz, ...
   ener,enef,enez,eber,ebef,ebez,eper,epef,epez]=magnetic(rho,theta);
else
  [b,dbds,dbdt,bp,sflx,dsdr,dsdz,dtdr,dtdz, ...
   ener,enef,enez,eber,ebef,ebez,eper,epef,epez,...
   dbds2,dbdst,dbdt2,dsdr2,dsdrz,dsdz2,dsdr3,dsdr2z,dsdrz2,dsdz3...
   dtdr2,dtdrz,dtdz2]=magnetic(rho,theta);
end
%
% ----- Wave vector and refraction index -------------------------------------
coom=cnst.c/plasma.omega; coomsq=coom^2; om2=plasma.omega^2;
kn=kr.*ener +kf.*enef +kz.*enez;  Nn=coom*kn;  Nn2=Nn.^2;
kb=kr.*eber +kf.*ebef +kz.*ebez;  Nb=coom*kb;  Nb2=Nb.^2;
kp=kr.*eper +kf.*epef +kz.*epez;  Np=coom*kp;  Np2=Np.^2; N2=Nn2+Nb2+Np2;
%
% ----- Parabolic profiles and logarithmic derivatives -----------------------
nspec=size(plasma.amass,2);
p=(1-sflx.^2*plasma.na);
for k=1:nspec
 dLNpds(:,k) =-2*sflx.*plasma.na(k).*plasma.nb(k)./p(:,k);
 dLNpds2(:,k)=dLNpds(:,k)./sflx-dLNpds(:,k).^2./plasma.nb(k);
 n(:,k)=plasma.n0(k).*p(:,k).^plasma.nb(k);                % Density
 T(:,k)=plasma.t0(k).*(1-sflx.^2*plasma.ta(k)).^plasma.nb(k);     % Temperature
end
dLNnds =dLNpds;
dLNnds2=dLNpds2;
%
% ----- Plasma parameters and logarithmic derivatives ------------------------
omp2= ones(size(theta))*((plasma.acharge'.*cnst.e).^2./...
    (plasma.amass'.*cnst.mp*cnst.eps0))'.*n;
omp = sqrt(omp2);                            % Plasma frequencies
omc = b*((plasma.acharge'.*cnst.e)./(plasma.amass'.*cnst.mp))';      % Cyclotron frequencies
omc2= omc.^2;
caoc2 = 1./sum(omp2./omc2, 2);               % (c_Alfven/c_light)^2
dLNomp2ds=dLNnds;                            % Logarithmic derivatives
one=ones(size(plasma.amass'))';
dLNomcds =(dbds./b)*one;
dLNomcdt =(dbdt./b)*one;
if (eval2nd)
 dLNomcds2=(dbds2./b-(dbds.^2   )./(b.^2))*one;
 dLNomcdst=(dbdst./b-(dbds.*dbdt)./(b.^2))*one;
 dLNomcdt2=(dbdt2./b-(dbdt.^2   )./(b.^2))*one;
 vth2=(2000.*T.*cnst.e)./(ones(size(T,1),1)*plasma.amass.*cnst.mp);
 vth=sqrt(vth2);                             % Thermal velocities
 nu=4.8E-14*15*n.*(1000*T).^(-1.5)./(ones(size(T,1),1)*sqrt(plasma.amass));% Collisions NRL p.34
 nu(1)=nu(1)*60.42;
end
%
switch plasma.MODEL(1:6)
%
% +++++ Cold plasma model ++++++++++++++++++++++++++++++++++++++++++++++++++++
%
case {'msw1x1','cld2x2','cld3x3'}

  % --- Elementary functions and derivatives ---------------------------------
  if (isCmplx && ~eval1st)
    iomceps=0.003*i*om2;
  else
    iomceps=0;
  end


  omc2Mom2= omc2-om2-2*iomceps;              % S,D,P
  Si=omp2./omc2Mom2;   S=1+sum( Si ,2);
  Di=(omc/plasma.omega).*Si;  D=  sum( Di ,2);
  Pi=omp2/om2;         P=1-sum( Pi ,2);

                                             % 1st derivatives S,D,P
  zroi=zeros(size(Si));                         zro=zeros(size(S));
  dLNSids  =dLNnds-2.*omc2./omc2Mom2.*dLNomcds; dSds = sum( Si.*dLNSids ,2);
  dLNSidt  =      -2.*omc2./omc2Mom2.*dLNomcdt; dSdt = sum( Si.*dLNSidt ,2);
  dLNSidom =2*plasma.omega./omc2Mom2;                  dSdom= sum( Si.*dLNSidom,2);
  dLNDids  =dLNSids + dLNomcds;                 dDds = sum( Di.*dLNDids ,2);
  dLNDidt  =dLNSidt + dLNomcdt;                 dDdt = sum( Di.*dLNDidt ,2);
  dLNDidom =(3*om2-omc2)./(plasma.omega*omc2Mom2);     dDdom= sum( Di.*dLNDidom,2);
  if (strcmp(plasma.MODEL(1:6),'cld3x3')||strcmp(oper,'Pol')||strcmp(oper,'Amp'))
   dLNPids =dLNnds;                             dPds =-sum( Pi.*dLNPids ,2);
   dLNPidt =zroi;                               dPdt = zro;
   dLNPidom=2/plasma.omega;                     dPdom= sum( Pi.*dLNPidom,2);
  end
                                             % 2nd derivatives S,D,P
  if (eval2nd)
   omOM=om2 ./omc2Mom2; 
   ocOM=omc2./omc2Mom2;
   dLNSids2= 2*ocOM.*(2*omOM.*dLNomcds.*dLNomcds -dLNomcds2)+dLNnds2;
   dLNSidst= 2*ocOM.*(2*omOM.*dLNomcds.*dLNomcdt -dLNomcdst);
   dLNSidt2= 2*ocOM.*(2*omOM.*dLNomcdt.*dLNomcdt -dLNomcdt2);
   dLNDids2= dLNomcds2 +dLNSids2;
   dLNDidst= dLNomcdst +dLNSidst;
   dLNDidt2= dLNomcdt2 +dLNSidt2;

   dSids2 = Si.*(dLNSids2 +dLNSids.*dLNSids);   dSds2=sum( dSids2 ,2);
   dSidst = Si.*(dLNSidst +dLNSids.*dLNSidt);   dSdst=sum( dSidst ,2);
   dSidt2 = Si.*(dLNSidt2 +dLNSidt.*dLNSidt);   dSdt2=sum( dSidt2 ,2);
   dDids2 = Di.*(dLNDids2 +dLNDids.*dLNDids);   dDds2=sum( dDids2 ,2);
   dDidst = Di.*(dLNDidst +dLNDids.*dLNDidt);   dDdst=sum( dDidst ,2);
   dDidt2 = Di.*(dLNDidt2 +dLNDidt.*dLNDidt);   dDdt2=sum( dDidt2 ,2);

   if (strcmp(plasma.MODEL(1:6),'cld3x3')||strcmp(oper,'Pol')||strcmp(oper,'Amp'))
    dLNPids2= dLNnds2;   dLNPidst= zroi;  dLNPidt2= zroi;
    dPids2= Pi.*(dLNPids2 +dLNPids.*dLNPids);   dPds2=sum( dPids2 ,2);
    dPidst= Pi.*(dLNPidst +dLNPids.*dLNPidt);   dPdst=sum( dPidst ,2);
    dPidt2= Pi.*(dLNPidt2 +dLNPidt.*dLNPidt);   dPdt2=sum( dPidt2 ,2);
   end;
  end;

  % --- Tensor elements and derivatives --------------------------------------

  D11 = Nb2+Np2-S;        dD11ds =-dSds;          dD11dt =-dSdt;
  D12 =-Nn.*Nb-D*i;       dD12ds =-dDds*i;        dD12dt =-dDdt*i;
  D22 = Nn2+Np2-S;        dD22ds =-dSds;          dD22dt =-dSdt;
  dD11dkn= zro;           dD11dkb= Nb*coom*2;     dD11dkp= Np*coom*2;
  dD12dkn=-Nb*coom;       dD12dkb=-Nn*coom;       dD12dkp= zro;
  dD22dkn= Nn*coom*2;     dD22dkb= zro;           dD22dkp= Np*coom*2;
  dD11dom=-2/plasma.omega*(Nb2+Np2)-dSdom;
  dD12dom=-2/plasma.omega*(Nn.*Nb) -dDdom*i;
  dD22dom=-2/plasma.omega*(Nn2+Np2)-dSdom;

  if (strcmp(plasma.MODEL(1:6),'cld3x3')||strcmp(oper,'Pol')||strcmp(oper,'Amp'))
   D13=-Nn.*Np;           dD13ds = zro;           dD13dt = zro;
   D23=-Nb.*Np;           dD23ds = zro;           dD23dt = zro;
   D33= Nn2+Nb2-P;        dD33ds =-dPds;          dD33dt =-dPdt; 
   dD13dkn=-Np*coom;      dD13dkb= zro;           dD13dkp=-Nn*coom;
   dD23dkn= zro;          dD23dkb=-Np*coom;       dD23dkp=-Nb*coom;
   dD33dkn= Nn*coom*2;    dD33dkb= Nb*coom*2;     dD33dkp= zro;
   dD13dom=-2/plasma.omega*(Nn.*Np);          dD23dom=-2/plasma.omega*(Nb.*Np);
   dD33dom=-2/plasma.omega*(Nn2+Nb2)-dPdom;
  end;

  if (eval2nd)
   dD11ds2=-dSds2;        dD11dst=-dSdst;         dD11dt2=-dSdt2;
   dD12ds2=-dDds2*i;      dD12dst=-dDdst*i;       dD12dt2=-dDdt2*i;
   dD22ds2=-dSds2;        dD22dst=-dSdst;         dD22dt2=-dSdt2;
   dD11dkn2 = zro;        dD11dkb2 = coomsq*2;    dD11dkp2 = coomsq*2;
   dD11dknkb= zro;        dD11dknkp= zro;         dD11dkbkp= zro;
   dD12dkn2 = zro;        dD12dkb2 = zro;         dD12dkp2 = zro;
   dD12dknkb=-coomsq;     dD12dknkp= zro;         dD12dkbkp= zro;
   dD22dkn2 = coomsq*2;   dD22dkb2 = zro;         dD22dkp2 = coomsq*2;
   dD22dknkb= zro;        dD22dknkp= zro;         dD22dkbkp= zro;

   if (strcmp(plasma.MODEL(1:6),'cld3x3'))
    dD13ds2 = zro;        dD13dst = zro;          dD13dt2 = zro;
    dD23ds2 = zro;        dD23dst = zro;          dD23dt2 = zro;
    dD33ds2 =-dPds2;      dD33dst =-dPdst;        dD33dt2 =-dPdt2;
    dD13dkn2 = zro;       dD13dkb2 = zro;         dD13dkp2 = zro;
    dD23dkn2 = zro;       dD23dkb2 = zro;         dD23dkp2 = zro;
    dD33dkn2 = coomsq*2;  dD33dkb2 = coomsq*2;    dD33dkp2 = zro;
    dD13dknkb= zro;       dD13dknkp=-coomsq;      dD13dkbkp= zro;
    dD23dknkb= zro;       dD23dknkp= zro;         dD23dkbkp=-coomsq;
    dD33dknkb= zro;       dD33dknkp= zro;         dD33dkbkp= zro;
   end;
  end

  % --- Additional quantities for symplectic matching ------------------------

  if (eval2nd)
   dD11dr =dD11ds.*dsdr + dD11dt.*dtdr; dD11f=zro;
   dD11dz =dD11ds.*dsdz + dD11dt.*dtdz;
   dD12dr =dD12ds.*dsdr + dD12dt.*dtdr; dD12f=zro;
   dD12dz =dD12ds.*dsdz + dD12dt.*dtdz;
   dD22dr =dD22ds.*dsdr + dD22dt.*dtdr; dD22f=zro;
   dD22dz =dD22ds.*dsdz + dD22dt.*dtdz; 

   dD11dkr=dD11dkn.*ener+dD11dkb.*eber+dD11dkp.*eper;
   dD11dkf=dD11dkn.*enef+dD11dkb.*ebef+dD11dkp.*epef;
   dD11dkz=dD11dkn.*enez+dD11dkb.*ebez+dD11dkp.*epez;
   dD12dkr=dD12dkn.*ener+dD12dkb.*eber+dD12dkp.*eper;
   dD12dkf=dD12dkn.*enef+dD12dkb.*ebef+dD12dkp.*epef;
   dD12dkz=dD12dkn.*enez+dD12dkb.*ebez+dD12dkp.*epez;
   dD22dkr=dD22dkn.*ener+dD22dkb.*eber+dD22dkp.*eper;
   dD22dkf=dD22dkn.*enef+dD22dkb.*ebef+dD22dkp.*epef;
   dD22dkz=dD22dkn.*enez+dD22dkb.*ebez+dD22dkp.*epez;

   gD11=cat(2,dD11dr,dD11dz,dD11dkr,dD11dkz);
   gD12=cat(2,dD12dr,dD12dz,dD12dkr,dD12dkz);
   gD22=cat(2,dD22dr,dD22dz,dD22dkr,dD22dkz);

                                             % Derivatives in osculating plane
   A= ydot(:,1).*yddot(:,3)+ydot(:,2).*yddot(:,4) ...
     -ydot(:,3).*yddot(:,1)-ydot(:,4).*yddot(:,2);
   sqmA=1./sqrt(abs(A));
   eq= sqmA*ones(1,4).* ydot(:,1:4);         % Velocity and acceleration
   ep= sqmA*ones(1,4).*yddot(:,1:4);
   dD11dq=sum(eq.*gD11,2); dD11dp=sum(ep.*gD11,2);
   dD12dq=sum(eq.*gD12,2); dD12dp=sum(ep.*gD12,2);
   dD22dq=sum(eq.*gD22,2); dD22dp=sum(ep.*gD22,2);

   if (strcmp(plasma.MODEL(1:6),'cld3x3'))
    dD13dr =dD13ds.* dsdr+dD13dt.* dtdr; dD13f =zro;
    dD13dz =dD13ds.* dsdz+dD13dt.* dtdz;
    dD23dr =dD23ds.* dsdr+dD23dt.* dtdr; dD23f =zro;
    dD23dz =dD23ds.* dsdz+dD23dt.* dtdz;
    dD33dr =dD33ds.* dsdr+dD33dt.* dtdr; dD33f =zro;
    dD33dz =dD33ds.* dsdz+dD33dt.* dtdz;
    dD13dkr=dD13dkn.*ener+dD13dkb.*eber+dD13dkp.*eper;
    dD13dkf=dD13dkn.*enef+dD13dkb.*ebef+dD13dkp.*epef;
    dD13dkz=dD13dkn.*enez+dD13dkb.*ebez+dD13dkp.*epez;
    dD23dkr=dD23dkn.*ener+dD23dkb.*eber+dD23dkp.*eper;
    dD23dkf=dD23dkn.*enef+dD23dkb.*ebef+dD23dkp.*epef;
    dD23dkz=dD23dkn.*enez+dD23dkb.*ebez+dD23dkp.*epez;
    dD33dkr=dD33dkn.*ener+dD33dkb.*eber+dD33dkp.*eper;
    dD33dkf=dD33dkn.*enef+dD33dkb.*ebef+dD33dkp.*epef;
    dD33dkz=dD33dkn.*enez+dD33dkb.*ebez+dD33dkp.*epez;
    gD13 =cat(2,dD13dr,dD13dz,dD13dkr,dD13dkz );
    gD23 =cat(2,dD23dr,dD23dz,dD23dkr,dD23dkz );
    gD33 =cat(2,dD33dr,dD33dz,dD33dkr,dD33dkz );
    dD13dq=sum(eq.*gD13,2); dD13dp=sum(ep.*gD13,2);
    dD23dq=sum(eq.*gD23,2); dD23dp=sum(ep.*gD23,2);
    dD33dq=sum(eq.*gD33,2); dD33dp=sum(ep.*gD33,2);
   end;
  end; % Con or Zro

% +++++ Warm plasma model ++++++++++++++++++++++++++++++++++++++++++++++++++++
otherwise
  error(['Unknown model ''' plasma.MODEL(1:6) '''.']);
end;
%
% ===== Dispersion function and derivatives ==================================
%
% +++++ Cold 1x1 fast magnetosonic wave approximation in 4D ++++++++++++++++++
if (strcmp(plasma.MODEL(1:6),'msw1x1') || strcmp(oper,'Msw'))
 U  =0.5*(1./caoc2-N2.*ones(size(theta)));   % Dispersion function
 trU = U;
 mon2 = caoc2.*U;                            % Monitor
 dUdom= 2*N2/plasma.omega;
 if (strcmp(oper,'Trj')||strcmp(oper,'Ant'))
  dUds=0.5*(2.*(dbds./b) -dLNpds(:,1))./caoc2;
  dUdt=0.5*(2.*(dbdt./b)             )./caoc2;
  dUdr = dUds.*dsdr + dUdt.*dtdr;
  dUdz = dUds.*dsdz + dUdt.*dtdz;
  dUdkn=coomsq * kn;
  dUdkb=coomsq * kb;
  dUdkp=coomsq * kp;
 end; %if(Trj,Ant)
end; % if(msw1x1)

% +++++ Cold 2x2 plasma dispersion relation ++++++++++++++++++++++++++++++++++

if ((strcmp(plasma.MODEL(1:6),'cld2x2')||strcmp(plasma.MODEL(1:6),'cld3x3')) && ...
    ~strcmp(oper,'Msw'))
 cD12=conj(D12);                             % Alias
 DD  =[D11 D12 cD12 D22];                    % Dispersion tensor
 U   = D11.*D22-D12.*cD12;                   % Determinant
 V   = 0.5*(D11+D22);        Vsq=V.^2;       % Trace/2
 T=U;  %mon2=V;                               % Trace is conversion monitor
 mon2 = log(abs(V)); % modification by Steve to see if it works better
 dm1o1=1./(Vsq-U);
 dp1o2=sqrt(Vsq-U);    dm1o2=1./dp1o2;
 dp3o2=dp1o2.*(Vsq-U); dm3o2=1./dp3o2;
 zro  =zeros(size(T));
% note: the sign pm is taken since along a ray U~0 -> dp1o2~V and we want
% the eigenvalue which is ~0.
 pm=-sign(V); eig2=V+pm.*dp1o2;              % Eigenvalue
 if (strcmp(oper,'Pol')||strcmp(oper,'Amp'))  % Eigenvector (2x2)
   dim2=sqrt(size(DD,2)); pol2=[]; pol3=[];
   for k=1:size(DD,1)                        %   for each ray
    dd=reshape(DD(k,:),dim2,dim2);
    new= ones(dim2,1); old=zeros(dim2,1); 
    n1=1; n2=1; tol=1E-4; itmax=100; it=0;
    while (n1>tol && n2>tol && it<=itmax)
      n1=norm((new-old)./new,inf);
      n2=norm((new+old)./new,inf);
      old=new; it=it+1; new=dd*new; new=new/norm(new);
    end
    if (it>itmax)
      rays.monErrAbort = 1; % hack to get my loop to continue
      error('dispertok:pol','Polarization not converged');
      %disp(['dispertok: polarization ' sprintf('ray#%i',k) ' not converged'])
      %old, new, pause
    end
    new3 = [new; -(D13(k)*new(1)+D23(k)*new(2))/D33(k)];
    pol2 = [pol2 new'];
    pol3 = [pol3 new3'];                     % Eigenvector (3x3 perturbation)
   end
   pol2=reshape(pol2',2,size(T,1))'; pol=pol2;
   pol3=reshape(pol3',3,size(T,1))';
 end

 dVdom=0.5*(dD11dom +dD22dom);
 dUdom=dD11dom.*D22 +D11.*dD22dom -2*real(dD12dom.*cD12);
 dTdom=dVdom +pm.*(V.*dVdom -0.5*dUdom).*dm1o2;

 if (eval1st)
   dUds =dD11ds .*D22 +D11.*dD22ds  -2*real(D12.*conj(dD12ds));   dVds =0.5*(dD11ds +dD22ds );
   dUdt =dD11dt .*D22 +D11.*dD22dt  -2*real(D12.*conj(dD12dt));   dVdt =0.5*(dD11dt +dD22dt );
   dUdkn=dD11dkn.*D22 +D11.*dD22dkn -2*real(D12.*conj(dD12dkn));  dVdkn=0.5*(dD11dkn+dD22dkn);
   dUdkb=dD11dkb.*D22 +D11.*dD22dkb -2*real(D12.*conj(dD12dkb));  dVdkb=0.5*(dD11dkb+dD22dkb);
   dUdkp=dD11dkp.*D22 +D11.*dD22dkp -2*real(D12.*conj(dD12dkp));  dVdkp=0.5*(dD11dkp+dD22dkp);
   % these derivatives are not used, rather dTdx, etc
%    dTds =dVds  +pm.*(V.*dVds  -0.5*dUds ).*dm1o2;
%    dTdt =dVdt  +pm.*(V.*dVdt  -0.5*dUdt ).*dm1o2;
%    dTdkn=dVdkn +pm.*(V.*dVdkn -0.5*dUdkn).*dm1o2;
%    dTdkb=dVdkb +pm.*(V.*dVdkb -0.5*dUdkb).*dm1o2;
%    dTdkp=dVdkp +pm.*(V.*dVdkp -0.5*dUdkp).*dm1o2;
 end; %if(Trj,Ant,Con,Zro,Amp)

 if (eval2nd)                                % Hessian and others
   H11 =2*dD11dq.*dD22dq                 -2*     dD12dq.*conj(dD12dq);
   H22 =2*dD11dp.*dD22dp                 -2*     dD12dp.*conj(dD12dp);
   H12 =  dD11dq.*dD22dp +dD11dp.*dD22dq -2*real(dD12dq.*conj(dD12dp));
   dUdq=  dD11dq.* D22    +D11  .*dD22dq -2*real( D12 .* conj(dD12dq));
   dUdp=  dD11dp.* D22    +D11  .*dD22dp -2*real( D12 .* conj(dD12dp));

   dUds2  =dD11ds2  .*D22+dD11ds .*dD22ds +dD11ds .*dD22ds +D11.*dD22ds2  -2*real(dD12ds .*conj(dD12ds) +D12.*conj(dD12ds2));
   dUdst  =dD11dst  .*D22+dD11ds .*dD22dt +dD11dt .*dD22ds +D11.*dD22dst  -2*real(dD12dt .*conj(dD12ds) +D12.*conj(dD12dst));
   dUdt2  =dD11dt2  .*D22+dD11dt .*dD22dt +dD11dt .*dD22dt +D11.*dD22dt2  -2*real(dD12dt .*conj(dD12dt) +D12.*conj(dD12dt2));
   dUdkn2 =dD11dkn2 .*D22+dD11dkn.*dD22dkn+dD11dkn.*dD22dkn+D11.*dD22dkn2 -2*real(dD12dkn.*conj(dD12dkn)+D12.*conj(dD12dkn2 ));
   dUdknkb=dD11dknkb.*D22+dD11dkn.*dD22dkb+dD11dkb.*dD22dkn+D11.*dD22dknkb-2*real(dD12dkb.*conj(dD12dkn)+D12.*conj(dD12dknkb));
   dUdknkp=dD11dknkp.*D22+dD11dkn.*dD22dkp+dD11dkp.*dD22dkn+D11.*dD22dknkp-2*real(dD12dkp.*conj(dD12dkn)+D12.*conj(dD12dknkp));
   dUdkb2 =dD11dkb2 .*D22+dD11dkb.*dD22dkb+dD11dkb.*dD22dkb+D11.*dD22dkb2 -2*real(dD12dkb.*conj(dD12dkb)+D12.*conj(dD12dkb2 ));
   dUdkbkp=dD11dkbkp.*D22+dD11dkb.*dD22dkp+dD11dkp.*dD22dkb+D11.*dD22dkbkp-2*real(dD12dkp.*conj(dD12dkb)+D12.*conj(dD12dkbkp));
   dUdkp2 =dD11dkp2 .*D22+dD11dkp.*dD22dkp+dD11dkp.*dD22dkp+D11.*dD22dkp2 -2*real(dD12dkp.*conj(dD12dkp)+D12.*conj(dD12dkp2 ));
   dUdskn =               dD11ds .*dD22dkn+dD11dkn.*dD22ds                -2*real(dD12dkn.*conj(dD12ds));
   dUdskb =               dD11ds .*dD22dkb+dD11dkb.*dD22ds                -2*real(dD12dkb.*conj(dD12ds));
   dUdskp =               dD11ds .*dD22dkp+dD11dkp.*dD22ds                -2*real(dD12dkp.*conj(dD12ds));
   dUdtkn =               dD11dt .*dD22dkn+dD11dkn.*dD22dt                -2*real(dD12dkn.*conj(dD12dt));
   dUdtkb =               dD11dt .*dD22dkb+dD11dkb.*dD22dt                -2*real(dD12dkb.*conj(dD12dt));
   dUdtkp =               dD11dt .*dD22dkp+dD11dkp.*dD22dt                -2*real(dD12dkp.*conj(dD12dt));

   dVds2  =0.5*(dD11ds2  +dD22ds2  );   dVdst  =0.5*(dD11dst  +dD22dst  );   dVdt2  =0.5*(dD11dt2  +dD22dt2  );
   dVdkn2 =0.5*(dD11dkn2 +dD22dkn2 );   dVdknkb=0.5*(dD11dknkb+dD22dknkb);   dVdknkp=0.5*(dD11dknkp+dD22dknkp);
   dVdkb2 =0.5*(dD11dkb2 +dD22dkb2 );   dVdkbkp=0.5*(dD11dkbkp+dD22dkbkp);
   dVdkp2 =0.5*(dD11dkp2 +dD22dkp2 );
   dVdskn =zro;                         dVdskb =zro;                         dVdskp =zro;                      
   dVdtkn =zro;                         dVdtkb =zro;                         dVdtkp =zro;                      

   % These derivs are not used, rather dTdx2, etc
%    dTds2  =dVds2   +pm.*dm1o2.*(-dm1o1.*(V.*dVds -0.5*dUds ).^2                    +dVds.^2      +V.*dVds2   -0.5*dUds2  );
%    dTdst  =dVdst   +pm.*dm1o2.*(-dm1o1.*(V.*dVds -0.5*dUds ).*(V.*dVdt -0.5*dUdt ) +dVds.* dVdt  +V.*dVdst   -0.5*dUdst  );
%    dTdt2  =dVdt2   +pm.*dm1o2.*(-dm1o1.*(V.*dVdt -0.5*dUdt ).^2                    +dVdt.^2      +V.*dVdt2   -0.5*dUdt2  );
%    dTdkn2 =dVdkn2  +pm.*dm1o2.*(-dm1o1.*(V.*dVdkn-0.5*dUdkn).^2                    +dVdkn.^2     +V.*dVdkn2  -0.5*dUdkn2 );
%    dTdknkb=dVdknkb +pm.*dm1o2.*(-dm1o1.*(V.*dVdkn-0.5*dUdkn).*(V.*dVdkb-0.5*dUdkb) +dVdkn.*dVdkb +V.*dVdknkb -0.5*dUdknkb);
%    dTdknkp=dVdknkp +pm.*dm1o2.*(-dm1o1.*(V.*dVdkn-0.5*dUdkn).*(V.*dVdkp-0.5*dUdkp) +dVdkn.*dVdkp +V.*dVdknkp -0.5*dUdknkp);
%    dTdkb2 =dVdkb2  +pm.*dm1o2.*(-dm1o1.*(V.*dVdkb-0.5*dUdkb).^2                    +dVdkb.^2     +V.*dVdkb2  -0.5*dUdkb2 );
%    dTdkbkp=dVdkbkp +pm.*dm1o2.*(-dm1o1.*(V.*dVdkb-0.5*dUdkb).*(V.*dVdkp-0.5*dUdkp) +dVdkb.*dVdkp +V.*dVdkbkp -0.5*dUdkbkp);
%    dTdkp2 =dVdkp2  +pm.*dm1o2.*(-dm1o1.*(V.*dVdkp-0.5*dUdkp).^2                    +dVdkp.^2     +V.*dVdkp2  -0.5*dUdkp2 );
%    dTdskn =dVdskn  +pm.*dm1o2.*(-dm1o1.*(V.*dVds -0.5*dUds ).*(V.*dVdkn-0.5*dUdkn) +dVds .*dVdkn +V.*dVdskn  -0.5*dUdskn );
%    dTdskb =dVdskb  +pm.*dm1o2.*(-dm1o1.*(V.*dVds -0.5*dUds ).*(V.*dVdkb-0.5*dUdkb) +dVds .*dVdkb +V.*dVdskb  -0.5*dUdskb );
%    dTdskp =dVdskp  +pm.*dm1o2.*(-dm1o1.*(V.*dVds -0.5*dUds ).*(V.*dVdkp-0.5*dUdkp) +dVds .*dVdkp +V.*dVdskp  -0.5*dUdskp );
%    dTdtkn =dVdtkn  +pm.*dm1o2.*(-dm1o1.*(V.*dVdt -0.5*dUdt ).*(V.*dVdkn-0.5*dUdkn) +dVdt .*dVdkn +V.*dVdtkn  -0.5*dUdtkn );
%    dTdtkb =dVdtkb  +pm.*dm1o2.*(-dm1o1.*(V.*dVdt -0.5*dUdt ).*(V.*dVdkb-0.5*dUdkb) +dVdt .*dVdkb +V.*dVdtkb  -0.5*dUdtkb );
%    dTdtkp =dVdtkp  +pm.*dm1o2.*(-dm1o1.*(V.*dVdt -0.5*dUdt ).*(V.*dVdkp-0.5*dUdkp) +dVdt .*dVdkp +V.*dVdtkp  -0.5*dUdtkp );

 end; % if(Con,Zro,Amp)

end; % if(cld2x2,cld3x3)

% +++++ Cold 3x3 plasma dispersion relation (extending 2x2) ++++++++++++++++++

if (strcmp(plasma.MODEL(1:6),'cld3x3') && ~strcmp(oper,'Msw'))
 U2      = U;      dU2dom  =dUdom;          % Alias for 2x2 quantities
 cD13 = conj(D13); cD23 = conj(D23); 
 DD=[D11 D12 D13 cD12 D22 D23 cD13 cD23 D33];% Dispersion tensor
                                            % Dispersion funtion
 %U = D33.*T2 -D23.*D11.*cD23 +D23.*D12.*cD13 -D13.*cD12.*cD23 +D13.*D22.*cD13;
 % (Steve) I think there is a sign error in the above expression...
 %    Also, T2 is U=det(2x2), I think.
 T2= D11.*D22-D12.*cD12;
 sgn_fix = -1;
 U = D33.*T2 -D23.*D11.*cD23 +D23.*D12.*cD13 ...
     +sgn_fix*(-D13.*cD12.*cD23 +D13.*D22.*cD13);
 
 MinrU = double(D22).*double(D33) -double(D23.*cD23) ...
        +double(D11).*double(D33) -double(D13.*cD13) ...
        +double(D11).*double(D22) -double(D12.*cD12);
 mon2= caoc2.*MinrU;
 dUdom= dD33dom.*U2 +D33.*dU2dom ...
   -dD23dom .*D11.*cD23 -D23     .*dD11dom .*cD23 -D23 .*D11.*conj(dD23dom) ...
   +dD23dom .*D12.*cD13 +D23     .*dD12dom .*cD13 +D23 .*D12.*conj(dD13dom) ...
   +sgn_fix*(-dD13dom.*cD12.*cD23 -D13.*conj(dD12dom).*cD23 -D13.*cD12.*conj(dD23dom)) ...
   +sgn_fix*(dD13dom .*D22.*cD13 +D13     .*dD22dom .*cD13 +D13 .*D22.*conj(dD13dom));

 if (eval1st)
  dU2ds   =dUds;    dU2dt   =dUdt;          % Alias for 2x2 quantities
  dU2dkn  =dUdkn;   dU2dkb  =dUdkb;   dU2dkp  =dUdkp;
  dUds = dD33ds .*T2 +D33.*dU2ds  ...
   -dD23ds  .*D11.*cD23 -D23     .*dD11ds  .*cD23 -D23 .*D11.*conj(dD23ds ) ...
   +dD23ds  .*D12.*cD13 +D23     .*dD12ds  .*cD13 +D23 .*D12.*conj(dD13ds ) ...
   +sgn_fix*(-dD13ds .*cD12.*cD23 -D13.*conj(dD12ds ).*cD23 -D13.*cD12.*conj(dD23ds )) ...
   +sgn_fix*(+dD13ds  .*D22.*cD13 +D13     .*dD22ds  .*cD13 +D13 .*D22.*conj(dD13ds ));
  dUdt = dD33dt .*T2 +D33.*dU2dt  ...
   -dD23dt  .*D11.*cD23 -D23     .*dD11dt  .*cD23 -D23 .*D11.*conj(dD23dt ) ...
   +dD23dt  .*D12.*cD13 +D23     .*dD12dt  .*cD13 +D23 .*D12.*conj(dD13dt ) ...
   +sgn_fix*(-dD13dt .*cD12.*cD23 -D13.*conj(dD12dt ).*cD23 -D13.*cD12.*conj(dD23dt )) ...
   +sgn_fix*(+dD13dt  .*D22.*cD13 +D13     .*dD22dt  .*cD13 +D13 .*D22.*conj(dD13dt ));
  dUdkn= dD33dkn.*T2 +D33.*dU2dkn ...
   -dD23dkn .*D11.*cD23 -D23     .*dD11dkn .*cD23 -D23 .*D11.*conj(dD23dkn) ...
   +dD23dkn .*D12.*cD13 +D23     .*dD12dkn .*cD13 +D23 .*D12.*conj(dD13dkn) ...
   +sgn_fix*(-dD13dkn.*cD12.*cD23 -D13.*conj(dD12dkn).*cD23 -D13.*cD12.*conj(dD23dkn)) ...
   +sgn_fix*(+dD13dkn .*D22.*cD13 +D13     .*dD22dkn .*cD13 +D13 .*D22.*conj(dD13dkn));
  dUdkb= dD33dkb.*T2 +D33.*dU2dkb ...
   -dD23dkb .*D11.*cD23 -D23     .*dD11dkb .*cD23 -D23 .*D11.*conj(dD23dkb) ...
   +dD23dkb .*D12.*cD13 +D23     .*dD12dkb .*cD13 +D23 .*D12.*conj(dD13dkb) ...
   +sgn_fix*(-dD13dkb.*cD12.*cD23 -D13.*conj(dD12dkb).*cD23 -D13.*cD12.*conj(dD23dkb)) ...
   +sgn_fix*(+dD13dkb .*D22.*cD13 +D13     .*dD22dkb .*cD13 +D13 .*D22.*conj(dD13dkb));
  dUdkp= dD33dkp.*T2 +D33.*dU2dkp ...
   -dD23dkp .*D11.*cD23 -D23     .*dD11dkp .*cD23 -D23 .*D11.*conj(dD23dkp) ...
   +dD23dkp .*D12.*cD13 +D23     .*dD12dkp .*cD13 +D23 .*D12.*conj(dD13dkp) ...
   +sgn_fix*(-dD13dkp.*cD12.*cD23 -D13.*conj(dD12dkp).*cD23 -D13.*cD12.*conj(dD23dkp)) ...
   +sgn_fix*(+dD13dkp .*D22.*cD13 +D13     .*dD22dkp .*cD13 +D13 .*D22.*conj(dD13dkp));
 end; %if(Trj,Ant,Con,Zro,Amp)

 %%%  FIX THIS: sgn_fix NOT ADDED BELOW!
 %%%  Trj should be okay, but Con, Zro, Amp are not
 
 % Aliases for 2x2 quantities
 if (eval2nd)
   cdD12dq=conj(dD12dq); cdD12dp=conj(dD12dp);
   cdD13dq=conj(dD13dq); cdD13dp=conj(dD13dp);
   cdD23dq=conj(dD23dq); cdD23dp=conj(dD23dp);
   dU2dq=dUdq;    dU2dp=dUdp;

   dU2ds2  =dUds2;   dU2dst  =dUdst;   dU2dt2  =dUdt2;
   dU2dskn =dUdskn;  dU2dskb =dUdskb;  dU2dskp =dUdskp;
   dU2dtkn =dUdtkn;  dU2dtkb =dUdtkb;  dU2dtkp =dUdtkp;
   dU2dkn2 =dUdkn2;  dU2dkb2 =dUdkb2;  dU2dkp2 =dUdkp2;
   dU2dknkb=dUdknkb; dU2dknkp=dUdknkp; dU2dkbkp=dUdkbkp; 

  dUdq= dD33dq.*T2 +D33.*dU2dq ...
   -dD23dq .*D11.*cD23 -D23     .*dD11dq .*cD23 -D23 .*D11.*conj(dD23dq) ...
   +dD23dq .*D12.*cD13 +D23     .*dD12dq .*cD13 +D23 .*D12.*conj(dD13dq) ...
   -dD13dq.*cD12.*cD23 -D13.*conj(dD12dq).*cD23 -D13.*cD12.*conj(dD23dq) ...
   +dD13dq .*D22.*cD13 +D13     .*dD22dq .*cD13 +D13 .*D22.*conj(dD13dq);
  dUdp= dD33dp.*T2 +D33.*dU2dp ...
   -dD23dp .*D11.*cD23 -D23     .*dD11dp .*cD23 -D23 .*D11.*conj(dD23dp) ...
   +dD23dp .*D12.*cD13 +D23     .*dD12dp .*cD13 +D23 .*D12.*conj(dD13dp) ...
   -dD13dp.*cD12.*cD23 -D13.*conj(dD12dp).*cD23 -D13.*cD12.*conj(dD23dp) ...
   +dD13dp .*D22.*cD13 +D13     .*dD22dp .*cD13 +D13 .*D22.*conj(dD13dp);

  H11= 2*(dD11dq .*D22 .*dD33dq +dD11dq.*dD22dq .*D33 +D11.* dD22dq.*dD33dq)...
      -2*(dD12dq.*cD12 .*dD33dq +dD12dq.*cdD12dq.*D33 +D12.*cdD12dq.*dD33dq)...
      +2*(dD13dq.*cD13 .*dD22dq +dD13dq.*cdD13dq.*D22 +D13.*cdD13dq.*dD22dq)...
      -2*(dD23dq.*cD23 .*dD11dq +dD23dq.*cdD23dq.*D11 +D23.*cdD23dq.*dD11dq)...
  +4*imag(dD12dq.*cD13 .*dD23dq +dD12dq.*cdD13dq.*D23 +D12.*cdD13dq.*dD23dq);

  H12=   (dD11dp .*D22 .*dD33dq +dD11dp.*dD22dq .*D33 +D11.* dD22dp.*dD33dq...
         +dD11dq .*D22 .*dD33dp +dD11dq.*dD22dp .*D33 +D11.* dD22dq.*dD33dp)...
        -(dD12dp.*cD12 .*dD33dq +dD12dp.*cdD12dq.*D33 +D12.*cdD12dp.*dD33dq...
         +dD12dq.*cD12 .*dD33dp +dD12dq.*cdD12dp.*D33 +D12.*cdD12dq.*dD33dp)...
        +(dD13dp.*cD13 .*dD22dq +dD13dp.*cdD13dq.*D22 +D13.*cdD13dp.*dD22dq...
         +dD13dq.*cD13 .*dD22dp +dD13dq.*cdD13dp.*D22 +D13.*cdD13dq.*dD22dp)...
        -(dD23dp.*cD23 .*dD11dq +dD23dp.*cdD23dq.*D11 +D23.*cdD23dp.*dD11dq...
         +dD23dq.*cD23 .*dD11dp +dD23dq.*cdD23dp.*D11 +D23.*cdD23dq.*dD11dp)...
  +2*imag(dD12dp.*cD13 .*dD23dq +dD12dp.*cdD13dq.*D23 +D12.*cdD13dp.*dD23dq...
         +dD12dq.*cD13 .*dD23dp +dD12dq.*cdD13dp.*D23 +D12.*cdD13dq.*dD23dp);

  H22= 2*(dD11dp .*D22 .*dD33dp +dD11dp.*dD22dp .*D33 +D11.* dD22dp.*dD33dp)...
      -2*(dD12dp.*cD12 .*dD33dp +dD12dp.*cdD12dp.*D33 +D12.*cdD12dp.*dD33dp)...
      +2*(dD13dp.*cD13 .*dD22dp +dD13dp.*cdD13dp.*D22 +D13.*cdD13dp.*dD22dp)...
      -2*(dD23dp.*cD23 .*dD11dp +dD23dp.*cdD23dp.*D11 +D23.*cdD23dp.*dD11dp)...
  +4*imag(dD12dp.*cD13 .*dD23dp +dD12dp.*cdD13dp.*D23 +D12.*cdD13dp.*dD23dp);

  dUds2= dD33ds2.*T2 +dD33ds.*dU2ds +dD33ds.*dU2ds +D33.*dU2ds2 ...
   -dD23ds2     .*D11     .*cD23   -dD23ds     .*dD11ds      .*cD23   -dD23ds    .*D11   .*conj(dD23ds)...
   +dD23ds2     .*D12     .*cD13   +dD23ds     .*dD12ds      .*cD13   +dD23ds    .*D12   .*conj(dD13ds)... 
   -dD13ds2    .*cD12     .*cD23   -dD13ds.*conj(dD12ds)     .*cD23   -dD13ds   .*cD12   .*conj(dD23ds)... 
   +dD13ds2     .*D22     .*cD13   +dD13ds     .*dD22ds      .*cD13   +dD13ds    .*D22   .*conj(dD13ds)...
   -dD23ds     .*dD11ds   .*cD23    -D23       .*dD11ds2     .*cD23    -D23     .*dD11ds .*conj(dD23ds)...
   +dD23ds     .*dD12ds   .*cD13    +D23       .*dD12ds2     .*cD13    +D23     .*dD12ds .*conj(dD13ds)...
   -dD13ds.*conj(dD12ds)  .*cD23    -D13  .*conj(dD12ds2)    .*cD23    -D13.*conj(dD12ds).*conj(dD23ds)...
   +dD13ds     .*dD22ds   .*cD13    +D13       .*dD22ds2     .*cD13    +D13     .*dD22ds .*conj(dD13ds)...
   -dD23ds      .*D11.*conj(dD23ds) -D23       .*dD11ds .*conj(dD23ds) -D23      .*D11   .*conj(dD23ds2)...
   +dD23ds      .*D12.*conj(dD13ds) +D23       .*dD12ds .*conj(dD13ds) +D23      .*D12   .*conj(dD13ds2)...
   -dD13ds     .*cD12.*conj(dD23ds) -D13  .*conj(dD12ds).*conj(dD23ds) -D13     .*cD12   .*conj(dD23ds2)...
   +dD13ds      .*D22.*conj(dD13ds) +D13       .*dD22ds .*conj(dD13ds) +D13      .*D22   .*conj(dD13ds2);

  dUdst= dD33dst.*T2 +dD33ds.*dU2dt +dD33dt.*dU2ds +D33.*dU2dst ...
   -dD23dst     .*D11     .*cD23   -dD23ds     .*dD11dt      .*cD23   -dD23ds    .*D11   .*conj(dD23dt)...
   +dD23dst     .*D12     .*cD13   +dD23ds     .*dD12dt      .*cD13   +dD23ds    .*D12   .*conj(dD13dt)... 
   -dD13dst    .*cD12     .*cD23   -dD13ds.*conj(dD12dt)     .*cD23   -dD13ds   .*cD12   .*conj(dD23dt)... 
   +dD13dst     .*D22     .*cD13   +dD13ds     .*dD22dt      .*cD13   +dD13ds    .*D22   .*conj(dD13dt)...
   -dD23dt     .*dD11ds   .*cD23    -D23       .*dD11dst     .*cD23    -D23     .*dD11ds .*conj(dD23dt)...
   +dD23dt     .*dD12ds   .*cD13    +D23       .*dD12dst     .*cD13    +D23     .*dD12ds .*conj(dD13dt)...
   -dD13dt.*conj(dD12ds)  .*cD23    -D13  .*conj(dD12dst)    .*cD23    -D13.*conj(dD12ds).*conj(dD23dt)...
   +dD13dt     .*dD22ds   .*cD13    +D13       .*dD22dst     .*cD13    +D13     .*dD22ds .*conj(dD13dt)...
   -dD23dt      .*D11.*conj(dD23ds) -D23       .*dD11dt .*conj(dD23ds) -D23      .*D11   .*conj(dD23dst)...
   +dD23dt      .*D12.*conj(dD13ds) +D23       .*dD12dt .*conj(dD13ds) +D23      .*D12   .*conj(dD13dst)...
   -dD13dt     .*cD12.*conj(dD23ds) -D13  .*conj(dD12dt).*conj(dD23ds) -D13     .*cD12   .*conj(dD23dst)...
   +dD13dt      .*D22.*conj(dD13ds) +D13       .*dD22dt .*conj(dD13ds) +D13      .*D22   .*conj(dD13dst);

  dUdt2= dD33dt2.*T2 +dD33dt.*dU2dt +dD33dt.*dU2dt +D33.*dU2dt2 ...
   -dD23dt2     .*D11     .*cD23   -dD23dt     .*dD11dt      .*cD23   -dD23dt    .*D11   .*conj(dD23dt)...
   +dD23dt2     .*D12     .*cD13   +dD23dt     .*dD12dt      .*cD13   +dD23dt    .*D12   .*conj(dD13dt)... 
   -dD13dt2    .*cD12     .*cD23   -dD13dt.*conj(dD12dt)     .*cD23   -dD13dt   .*cD12   .*conj(dD23dt)... 
   +dD13dt2     .*D22     .*cD13   +dD13dt     .*dD22dt      .*cD13   +dD13dt    .*D22   .*conj(dD13dt)...
   -dD23dt     .*dD11dt   .*cD23    -D23       .*dD11dt2     .*cD23    -D23     .*dD11dt .*conj(dD23dt)...
   +dD23dt     .*dD12dt   .*cD13    +D23       .*dD12dt2     .*cD13    +D23     .*dD12dt .*conj(dD13dt)...
   -dD13dt.*conj(dD12dt)  .*cD23    -D13  .*conj(dD12dt2)    .*cD23    -D13.*conj(dD12dt).*conj(dD23dt)...
   +dD13dt     .*dD22dt   .*cD13    +D13       .*dD22dt2     .*cD13    +D13     .*dD22dt .*conj(dD13dt)...
   -dD23dt      .*D11.*conj(dD23dt) -D23       .*dD11dt .*conj(dD23dt) -D23      .*D11   .*conj(dD23dt2)...
   +dD23dt      .*D12.*conj(dD13dt) +D23       .*dD12dt .*conj(dD13dt) +D23      .*D12   .*conj(dD13dt2)...
   -dD13dt     .*cD12.*conj(dD23dt) -D13  .*conj(dD12dt).*conj(dD23dt) -D13     .*cD12   .*conj(dD23dt2)...
   +dD13dt      .*D22.*conj(dD13dt) +D13       .*dD22dt .*conj(dD13dt) +D13      .*D22   .*conj(dD13dt2);

  dUdskn= dD33ds.*dU2dkn +dD33dkn.*dU2ds +D33.*dU2dskn ...
                                    -dD23ds     .*dD11dkn      .*cD23   -dD23ds    .*D11   .*conj(dD23dkn)...
                                    +dD23ds     .*dD12dkn      .*cD13   +dD23ds    .*D12   .*conj(dD13dkn)... 
                                    -dD13ds.*conj(dD12dkn)     .*cD23   -dD13ds   .*cD12   .*conj(dD23dkn)... 
                                    +dD13ds     .*dD22dkn      .*cD13   +dD13ds    .*D22   .*conj(dD13dkn)...
   -dD23dkn     .*dD11ds   .*cD23                                        -D23     .*dD11ds .*conj(dD23dkn)...
   +dD23dkn     .*dD12ds   .*cD13                                        +D23     .*dD12ds .*conj(dD13dkn)...
   -dD13dkn.*conj(dD12ds)  .*cD23                                        -D13.*conj(dD12ds).*conj(dD23dkn)...
   +dD13dkn     .*dD22ds   .*cD13                                        +D13     .*dD22ds .*conj(dD13dkn)...
   -dD23dkn      .*D11.*conj(dD23ds) -D23       .*dD11dkn .*conj(dD23ds)...
   +dD23dkn      .*D12.*conj(dD13ds) +D23       .*dD12dkn .*conj(dD13ds)...
   -dD13dkn     .*cD12.*conj(dD23ds) -D13  .*conj(dD12dkn).*conj(dD23ds)...
   +dD13dkn      .*D22.*conj(dD13ds) +D13       .*dD22dkn .*conj(dD13ds);

  dUdskb= dD33ds.*dU2dkb +dD33dkb.*dU2ds +D33.*dU2dskb ...
                                    -dD23ds     .*dD11dkb      .*cD23   -dD23ds    .*D11   .*conj(dD23dkb)...
                                    +dD23ds     .*dD12dkb      .*cD13   +dD23ds    .*D12   .*conj(dD13dkb)... 
                                    -dD13ds.*conj(dD12dkb)     .*cD23   -dD13ds   .*cD12   .*conj(dD23dkb)... 
                                    +dD13ds     .*dD22dkb      .*cD13   +dD13ds    .*D22   .*conj(dD13dkb)...
   -dD23dkb     .*dD11ds   .*cD23                                        -D23     .*dD11ds .*conj(dD23dkb)...
   +dD23dkb     .*dD12ds   .*cD13                                        +D23     .*dD12ds .*conj(dD13dkb)...
   -dD13dkb.*conj(dD12ds)  .*cD23                                        -D13.*conj(dD12ds).*conj(dD23dkb)...
   +dD13dkb     .*dD22ds   .*cD13                                        +D13     .*dD22ds .*conj(dD13dkb)...
   -dD23dkb      .*D11.*conj(dD23ds) -D23       .*dD11dkb .*conj(dD23ds)...
   +dD23dkb      .*D12.*conj(dD13ds) +D23       .*dD12dkb .*conj(dD13ds)...
   -dD13dkb     .*cD12.*conj(dD23ds) -D13  .*conj(dD12dkb).*conj(dD23ds)...
   +dD13dkb      .*D22.*conj(dD13ds) +D13       .*dD22dkb .*conj(dD13ds);

  dUdskp= dD33ds.*dU2dkp +dD33dkp.*dU2ds +D33.*dU2dskp ...
                                    -dD23ds     .*dD11dkp      .*cD23   -dD23ds    .*D11   .*conj(dD23dkp)...
                                    +dD23ds     .*dD12dkp      .*cD13   +dD23ds    .*D12   .*conj(dD13dkp)... 
                                    -dD13ds.*conj(dD12dkp)     .*cD23   -dD13ds   .*cD12   .*conj(dD23dkp)... 
                                    +dD13ds     .*dD22dkp      .*cD13   +dD13ds    .*D22   .*conj(dD13dkp)...
   -dD23dkp     .*dD11ds   .*cD23                                        -D23     .*dD11ds .*conj(dD23dkp)...
   +dD23dkp     .*dD12ds   .*cD13                                        +D23     .*dD12ds .*conj(dD13dkp)...
   -dD13dkp.*conj(dD12ds)  .*cD23                                        -D13.*conj(dD12ds).*conj(dD23dkp)...
   +dD13dkp     .*dD22ds   .*cD13                                        +D13     .*dD22ds .*conj(dD13dkp)...
   -dD23dkp      .*D11.*conj(dD23ds) -D23       .*dD11dkp .*conj(dD23ds)...
   +dD23dkp      .*D12.*conj(dD13ds) +D23       .*dD12dkp .*conj(dD13ds)...
   -dD13dkp     .*cD12.*conj(dD23ds) -D13  .*conj(dD12dkp).*conj(dD23ds)...
   +dD13dkp      .*D22.*conj(dD13ds) +D13       .*dD22dkp .*conj(dD13ds);

  dUdtkn= dD33dt.*dU2dkn +dD33dkn.*dU2dt +D33.*dU2dtkn ...
                                    -dD23dt     .*dD11dkn      .*cD23   -dD23dt    .*D11   .*conj(dD23dkn)...
                                    +dD23dt     .*dD12dkn      .*cD13   +dD23dt    .*D12   .*conj(dD13dkn)... 
                                    -dD13dt.*conj(dD12dkn)     .*cD23   -dD13dt   .*cD12   .*conj(dD23dkn)... 
                                    +dD13dt     .*dD22dkn      .*cD13   +dD13dt    .*D22   .*conj(dD13dkn)...
   -dD23dkn     .*dD11dt   .*cD23                                        -D23     .*dD11dt .*conj(dD23dkn)...
   +dD23dkn     .*dD12dt   .*cD13                                        +D23     .*dD12dt .*conj(dD13dkn)...
   -dD13dkn.*conj(dD12dt)  .*cD23                                        -D13.*conj(dD12dt).*conj(dD23dkn)...
   +dD13dkn     .*dD22dt   .*cD13                                        +D13     .*dD22dt .*conj(dD13dkn)...
   -dD23dkn      .*D11.*conj(dD23dt) -D23       .*dD11dkn .*conj(dD23dt)...
   +dD23dkn      .*D12.*conj(dD13dt) +D23       .*dD12dkn .*conj(dD13dt)...
   -dD13dkn     .*cD12.*conj(dD23dt) -D13  .*conj(dD12dkn).*conj(dD23dt)...
   +dD13dkn      .*D22.*conj(dD13dt) +D13       .*dD22dkn .*conj(dD13dt);

  dUdtkb= dD33dt.*dU2dkb +dD33dkb.*dU2dt +D33.*dU2dtkb ...
                                    -dD23dt     .*dD11dkb      .*cD23   -dD23dt    .*D11   .*conj(dD23dkb)...
                                    +dD23dt     .*dD12dkb      .*cD13   +dD23dt    .*D12   .*conj(dD13dkb)... 
                                    -dD13dt.*conj(dD12dkb)     .*cD23   -dD13dt   .*cD12   .*conj(dD23dkb)... 
                                    +dD13dt     .*dD22dkb      .*cD13   +dD13dt    .*D22   .*conj(dD13dkb)...
   -dD23dkb     .*dD11dt   .*cD23                                        -D23     .*dD11dt .*conj(dD23dkb)...
   +dD23dkb     .*dD12dt   .*cD13                                        +D23     .*dD12dt .*conj(dD13dkb)...
   -dD13dkb.*conj(dD12dt)  .*cD23                                        -D13.*conj(dD12dt).*conj(dD23dkb)...
   +dD13dkb     .*dD22dt   .*cD13                                        +D13     .*dD22dt .*conj(dD13dkb)...
   -dD23dkb      .*D11.*conj(dD23dt) -D23       .*dD11dkb .*conj(dD23dt)...
   +dD23dkb      .*D12.*conj(dD13dt) +D23       .*dD12dkb .*conj(dD13dt)...
   -dD13dkb     .*cD12.*conj(dD23dt) -D13  .*conj(dD12dkb).*conj(dD23dt)...
   +dD13dkb      .*D22.*conj(dD13dt) +D13       .*dD22dkb .*conj(dD13dt);

  dUdtkp= dD33dt.*dU2dkp +dD33dkp.*dU2dt +D33.*dU2dtkp ...
                                    -dD23dt     .*dD11dkp      .*cD23   -dD23dt    .*D11   .*conj(dD23dkp)...
                                    +dD23dt     .*dD12dkp      .*cD13   +dD23dt    .*D12   .*conj(dD13dkp)... 
                                    -dD13dt.*conj(dD12dkp)     .*cD23   -dD13dt   .*cD12   .*conj(dD23dkp)... 
                                    +dD13dt     .*dD22dkp      .*cD13   +dD13dt    .*D22   .*conj(dD13dkp)...
   -dD23dkp     .*dD11dt   .*cD23                                        -D23     .*dD11dt .*conj(dD23dkp)...
   +dD23dkp     .*dD12dt   .*cD13                                        +D23     .*dD12dt .*conj(dD13dkp)...
   -dD13dkp.*conj(dD12dt)  .*cD23                                        -D13.*conj(dD12dt).*conj(dD23dkp)...
   +dD13dkp     .*dD22dt   .*cD13                                        +D13     .*dD22dt .*conj(dD13dkp)...
   -dD23dkp      .*D11.*conj(dD23dt) -D23       .*dD11dkp .*conj(dD23dt)...
   +dD23dkp      .*D12.*conj(dD13dt) +D23       .*dD12dkp .*conj(dD13dt)...
   -dD13dkp     .*cD12.*conj(dD23dt) -D13  .*conj(dD12dkp).*conj(dD23dt)...
   +dD13dkp      .*D22.*conj(dD13dt) +D13       .*dD22dkp .*conj(dD13dt);

  dUdkn2= dD33dkn2.*T2 +dD33dkn.*dU2dkn +dD33dkn.*dU2dkn +D33.*dU2dkn2 ...
   -dD23dkn2     .*D11     .*cD23    -dD23dkn     .*dD11dkn      .*cD23    -dD23dkn   .*D11    .*conj(dD23dkn)...
   +dD23dkn2     .*D12.*     cD13    +dD23dkn     .*dD12dkn      .*cD13    +dD23dkn   .*D12    .*conj(dD13dkn)... 
   -dD13dkn2    .*cD12     .*cD23    -dD13dkn.*conj(dD12dkn)     .*cD23    -dD13dkn  .*cD12    .*conj(dD23dkn)... 
   +dD13dkn2     .*D22     .*cD13    +dD13dkn     .*dD22dkn      .*cD13    +dD13dkn   .*D22    .*conj(dD13dkn)...
   -dD23dkn     .*dD11dkn  .*cD23     -D23        .*dD11dkn2     .*cD23     -D23     .*dD11dkn .*conj(dD23dkn)...
   +dD23dkn     .*dD12dkn  .*cD13     +D23        .*dD12dkn2     .*cD13     +D23     .*dD12dkn .*conj(dD13dkn)...
   -dD13dkn.*conj(dD12dkn) .*cD23     -D13   .*conj(dD12dkn2)    .*cD23     -D13.*conj(dD12dkn).*conj(dD23dkn)...
   +dD13dkn     .*dD22dkn  .*cD13     +D13        .*dD22dkn2     .*cD13     +D13     .*dD22dkn .*conj(dD13dkn)...
   -dD23dkn      .*D11.*conj(dD23dkn) -D23        .*dD11dkn .*conj(dD23dkn) -D23      .*D11    .*conj(dD23dkn2)...
   +dD23dkn      .*D12.*conj(dD13dkn) +D23        .*dD12dkn .*conj(dD13dkn) +D23      .*D12    .*conj(dD13dkn2)...
   -dD13dkn     .*cD12.*conj(dD23dkn) -D13   .*conj(dD12dkn).*conj(dD23dkn) -D13     .*cD12    .*conj(dD23dkn2)...
   +dD13dkn      .*D22.*conj(dD13dkn) +D13        .*dD22dkn .*conj(dD13dkn) +D13      .*D22    .*conj(dD13dkn2);

  dUdkb2= dD33dkb2.*T2 +dD33dkb.*dU2dkb +dD33dkb.*dU2dkb +D33.*dU2dkb2 ...
   -dD23dkb2     .*D11     .*cD23    -dD23dkb     .*dD11dkb      .*cD23    -dD23dkb   .*D11    .*conj(dD23dkb)...
   +dD23dkb2     .*D12.*     cD13    +dD23dkb     .*dD12dkb      .*cD13    +dD23dkb   .*D12    .*conj(dD13dkb)... 
   -dD13dkb2    .*cD12     .*cD23    -dD13dkb.*conj(dD12dkb)     .*cD23    -dD13dkb  .*cD12    .*conj(dD23dkb)... 
   +dD13dkb2     .*D22     .*cD13    +dD13dkb     .*dD22dkb      .*cD13    +dD13dkb   .*D22    .*conj(dD13dkb)...
   -dD23dkb     .*dD11dkb  .*cD23     -D23        .*dD11dkb2     .*cD23     -D23     .*dD11dkb .*conj(dD23dkb)...
   +dD23dkb     .*dD12dkb  .*cD13     +D23        .*dD12dkb2     .*cD13     +D23     .*dD12dkb .*conj(dD13dkb)...
   -dD13dkb.*conj(dD12dkb) .*cD23     -D13   .*conj(dD12dkb2)    .*cD23     -D13.*conj(dD12dkb).*conj(dD23dkb)...
   +dD13dkb     .*dD22dkb  .*cD13     +D13        .*dD22dkb2     .*cD13     +D13     .*dD22dkb .*conj(dD13dkb)...
   -dD23dkb      .*D11.*conj(dD23dkb) -D23        .*dD11dkb .*conj(dD23dkb) -D23      .*D11    .*conj(dD23dkb2)...
   +dD23dkb      .*D12.*conj(dD13dkb) +D23        .*dD12dkb .*conj(dD13dkb) +D23      .*D12    .*conj(dD13dkb2)...
   -dD13dkb     .*cD12.*conj(dD23dkb) -D13   .*conj(dD12dkb).*conj(dD23dkb) -D13     .*cD12    .*conj(dD23dkb2)...
   +dD13dkb      .*D22.*conj(dD13dkb) +D13        .*dD22dkb .*conj(dD13dkb) +D13      .*D22    .*conj(dD13dkb2);

  dUdkp2= dD33dkp2.*T2 +dD33dkp.*dU2dkp +dD33dkp.*dU2dkp +D33.*dU2dkp2 ...
   -dD23dkp2     .*D11     .*cD23    -dD23dkp     .*dD11dkp      .*cD23    -dD23dkp   .*D11    .*conj(dD23dkp)...
   +dD23dkp2     .*D12.*     cD13    +dD23dkp     .*dD12dkp      .*cD13    +dD23dkp   .*D12    .*conj(dD13dkp)... 
   -dD13dkp2    .*cD12     .*cD23    -dD13dkp.*conj(dD12dkp)     .*cD23    -dD13dkp  .*cD12    .*conj(dD23dkp)... 
   +dD13dkp2     .*D22     .*cD13    +dD13dkp     .*dD22dkp      .*cD13    +dD13dkp   .*D22    .*conj(dD13dkp)...
   -dD23dkp     .*dD11dkp  .*cD23     -D23        .*dD11dkp2     .*cD23     -D23     .*dD11dkp .*conj(dD23dkp)...
   +dD23dkp     .*dD12dkp  .*cD13     +D23        .*dD12dkp2     .*cD13     +D23     .*dD12dkp .*conj(dD13dkp)...
   -dD13dkp.*conj(dD12dkp) .*cD23     -D13   .*conj(dD12dkp2)    .*cD23     -D13.*conj(dD12dkp).*conj(dD23dkp)...
   +dD13dkp     .*dD22dkp  .*cD13     +D13        .*dD22dkp2     .*cD13     +D13     .*dD22dkp .*conj(dD13dkp)...
   -dD23dkp      .*D11.*conj(dD23dkp) -D23        .*dD11dkp .*conj(dD23dkp) -D23      .*D11    .*conj(dD23dkp2)...
   +dD23dkp      .*D12.*conj(dD13dkp) +D23        .*dD12dkp .*conj(dD13dkp) +D23      .*D12    .*conj(dD13dkp2)...
   -dD13dkp     .*cD12.*conj(dD23dkp) -D13   .*conj(dD12dkp).*conj(dD23dkp) -D13     .*cD12    .*conj(dD23dkp2)...
   +dD13dkp      .*D22.*conj(dD13dkp) +D13        .*dD22dkp .*conj(dD13dkp) +D13      .*D22    .*conj(dD13dkp2);

  dUdknkb= dD33dknkb.*T2 +dD33dkn.*dU2dkb +dD33dkb.*dU2dkn +D33.*dU2dknkb ...
   -dD23dknkb    .*D11     .*cD23    -dD23dkn     .*dD11dkb      .*cD23    -dD23dkn   .*D11    .*conj(dD23dkb)...
   +dD23dknkb    .*D12.*     cD13    +dD23dkn     .*dD12dkb      .*cD13    +dD23dkn   .*D12    .*conj(dD13dkb)... 
   -dD13dknkb   .*cD12     .*cD23    -dD13dkn.*conj(dD12dkb)     .*cD23    -dD13dkn  .*cD12    .*conj(dD23dkb)... 
   +dD13dknkb    .*D22     .*cD13    +dD13dkn     .*dD22dkb      .*cD13    +dD13dkn   .*D22    .*conj(dD13dkb)...
   -dD23dkb     .*dD11dkn  .*cD23     -D23        .*dD11dknkb    .*cD23     -D23     .*dD11dkn .*conj(dD23dkb)...
   +dD23dkb     .*dD12dkn  .*cD13     +D23        .*dD12dknkb    .*cD13     +D23     .*dD12dkn .*conj(dD13dkb)...
   -dD13dkb.*conj(dD12dkn) .*cD23     -D13   .*conj(dD12dknkb)   .*cD23     -D13.*conj(dD12dkn).*conj(dD23dkb)...
   +dD13dkb     .*dD22dkn  .*cD13     +D13        .*dD22dknkb    .*cD13     +D13     .*dD22dkn .*conj(dD13dkb)...
   -dD23dkb      .*D11.*conj(dD23dkn) -D23        .*dD11dkb .*conj(dD23dkn) -D23      .*D11    .*conj(dD23dknkb)...
   +dD23dkb      .*D12.*conj(dD13dkn) +D23        .*dD12dkb .*conj(dD13dkn) +D23      .*D12    .*conj(dD13dknkb)...
   -dD13dkb     .*cD12.*conj(dD23dkn) -D13   .*conj(dD12dkb).*conj(dD23dkn) -D13     .*cD12    .*conj(dD23dknkb) ...
   +dD13dkb      .*D22.*conj(dD13dkn) +D13        .*dD22dkb .*conj(dD13dkn) +D13      .*D22    .*conj(dD13dknkb);

  dUdknkp= dD33dknkp.*T2 +dD33dkn.*dU2dkp +dD33dkp.*dU2dkn +D33.*dU2dknkp ...
   -dD23dknkp    .*D11     .*cD23    -dD23dkn     .*dD11dkp      .*cD23    -dD23dkn   .*D11    .*conj(dD23dkp)...
   +dD23dknkp    .*D12.*     cD13    +dD23dkn     .*dD12dkp      .*cD13    +dD23dkn   .*D12    .*conj(dD13dkp)... 
   -dD13dknkp   .*cD12     .*cD23    -dD13dkn.*conj(dD12dkp)     .*cD23    -dD13dkn  .*cD12    .*conj(dD23dkp)... 
   +dD13dknkp    .*D22     .*cD13    +dD13dkn     .*dD22dkp      .*cD13    +dD13dkn   .*D22    .*conj(dD13dkp)...
   -dD23dkp     .*dD11dkn  .*cD23     -D23        .*dD11dknkp    .*cD23     -D23     .*dD11dkn .*conj(dD23dkp)...
   +dD23dkp     .*dD12dkn  .*cD13     +D23        .*dD12dknkp    .*cD13     +D23     .*dD12dkn .*conj(dD13dkp)...
   -dD13dkp.*conj(dD12dkn) .*cD23     -D13   .*conj(dD12dknkp)   .*cD23     -D13.*conj(dD12dkn).*conj(dD23dkp)...
   +dD13dkp     .*dD22dkn  .*cD13     +D13        .*dD22dknkp    .*cD13     +D13     .*dD22dkn .*conj(dD13dkp)...
   -dD23dkp      .*D11.*conj(dD23dkn) -D23        .*dD11dkp .*conj(dD23dkn) -D23      .*D11    .*conj(dD23dknkp)...
   +dD23dkp      .*D12.*conj(dD13dkn) +D23        .*dD12dkp .*conj(dD13dkn) +D23      .*D12    .*conj(dD13dknkp)...
   -dD13dkp     .*cD12.*conj(dD23dkn) -D13   .*conj(dD12dkp).*conj(dD23dkn) -D13     .*cD12    .*conj(dD23dknkp)...
   +dD13dkp      .*D22.*conj(dD13dkn) +D13        .*dD22dkp .*conj(dD13dkn) +D13      .*D22    .*conj(dD13dknkp);

  dUdkbkp= dD33dkbkp.*T2 +dD33dkb.*dU2dkp +dD33dkp.*dU2dkb +D33.*dU2dkbkp ...
   -dD23dkbkp    .*D11     .*cD23    -dD23dkb     .*dD11dkp      .*cD23    -dD23dkb   .*D11    .*conj(dD23dkp)...
   +dD23dkbkp    .*D12.*     cD13    +dD23dkb     .*dD12dkp      .*cD13    +dD23dkb   .*D12    .*conj(dD13dkp)... 
   -dD13dkbkp   .*cD12     .*cD23    -dD13dkb.*conj(dD12dkp)     .*cD23    -dD13dkb  .*cD12    .*conj(dD23dkp)... 
   +dD13dkbkp    .*D22     .*cD13    +dD13dkb     .*dD22dkp      .*cD13    +dD13dkb   .*D22    .*conj(dD13dkp)...
   -dD23dkp     .*dD11dkb  .*cD23     -D23        .*dD11dkbkp    .*cD23     -D23     .*dD11dkb .*conj(dD23dkp)...
   +dD23dkp     .*dD12dkb  .*cD13     +D23        .*dD12dkbkp    .*cD13     +D23     .*dD12dkb .*conj(dD13dkp)...
   -dD13dkp.*conj(dD12dkb) .*cD23     -D13   .*conj(dD12dkbkp)   .*cD23     -D13.*conj(dD12dkb).*conj(dD23dkp)...
   +dD13dkp     .*dD22dkb  .*cD13     +D13        .*dD22dkbkp    .*cD13     +D13     .*dD22dkb .*conj(dD13dkp)...
   -dD23dkp      .*D11.*conj(dD23dkb) -D23        .*dD11dkp .*conj(dD23dkb) -D23      .*D11    .*conj(dD23dkbkp)...
   +dD23dkp      .*D12.*conj(dD13dkb) +D23        .*dD12dkp .*conj(dD13dkb) +D23      .*D12    .*conj(dD13dkbkp)...
   -dD13dkp     .*cD12.*conj(dD23dkb) -D13   .*conj(dD12dkp).*conj(dD23dkb) -D13     .*cD12    .*conj(dD23dkbkp)...
   +dD13dkp      .*D22.*conj(dD13dkb) +D13        .*dD22dkp .*conj(dD13dkb) +D13      .*D22    .*conj(dD13dkbkp);
 end; %if(Con,Zro,Amp)

end; %if(cld3x3)

% ----- Transform derivatives back to (r,f,z) ---------------------------------

if (eval1st)
 dUdr =dUds.*dsdr + dUdt.*dtdr;                        dVdr =dVds.*dsdr + dVdt.*dtdr;
 dUdz =dUds.*dsdz + dUdt.*dtdz;                        dVdz =dVds.*dsdz + dVdt.*dtdz;
 dUdkr=dUdkn.*ener +dUdkb.*eber +dUdkp.*eper;          dVdkr=dVdkn.*ener +dVdkb.*eber +dVdkp.*eper;
 dUdkf=dUdkn.*enef +dUdkb.*ebef +dUdkp.*epef;          dVdkf=dVdkn.*enef +dVdkb.*ebef +dVdkp.*epef;
 dUdkz=dUdkn.*enez +dUdkb.*ebez +dUdkp.*epez;          dVdkz=dVdkn.*enez +dVdkb.*ebez +dVdkp.*epez;
 dTdr =dVdr +pm.*(V.*dVdr -0.5*dUdr ).*dm1o2;
 dTdz =dVdz +pm.*(V.*dVdz -0.5*dUdz ).*dm1o2;
 dTdkr=dVdkr+pm.*(V.*dVdkr-0.5*dUdkr).*dm1o2;
 dTdkf=dVdkf+pm.*(V.*dVdkf-0.5*dUdkf).*dm1o2;
 dTdkz=dVdkz+pm.*(V.*dVdkz-0.5*dUdkz).*dm1o2;
end

if (eval2nd)                                 % --- Conversion + amplitude terms
 dUdrs=dUds2.*dsdr +dUdst.*dtdr;                       dVdrs=dVds2.*dsdr +dVdst.*dtdr;
 dUdzs=dUds2.*dsdz +dUdst.*dtdz;                       dVdzs=dVds2.*dsdz +dVdst.*dtdz;
 dUdrt=dUdst.*dsdr +dUdt2.*dtdr;                       dVdrt=dVdst.*dsdr +dVdt2.*dtdr;
 dUdzt=dUdst.*dsdz +dUdt2.*dtdz;                       dVdzt=dVdst.*dsdz +dVdt2.*dtdz;
 dUdskr=dUdskn  .*ener +dUdskb .*eber +dUdskp .*eper;  dVdskr=dVdskn  .*ener +dVdskb .*eber +dVdskp .*eper;
 %dUdskf=dUdskn  .*enef +dUdskb .*ebef +dUdskp .*epef;  dVdskf=dVdskn  .*enef +dVdskb .*ebef +dVdskp .*epef;
 dUdskz=dUdskn  .*enez +dUdskb .*ebez +dUdskp .*epez;  dVdskz=dVdskn  .*enez +dVdskb .*ebez +dVdskp .*epez;
 dUdtkr=dUdtkn  .*ener +dUdtkb .*eber +dUdtkp .*eper;  dVdtkr=dVdtkn  .*ener +dVdtkb .*eber +dVdtkp .*eper;
 %dUdtkf=dUdtkn  .*enef +dUdtkb .*ebef +dUdtkp .*epef;  dVdtkf=dVdtkn  .*enef +dVdtkb .*ebef +dVdtkp .*epef;
 dUdtkz=dUdtkn  .*enez +dUdtkb .*ebez +dUdtkp .*epez;  dVdtkz=dVdtkn  .*enez +dVdtkb .*ebez +dVdtkp .*epez;
 dUdkrkn=dUdkn2 .*ener +dUdknkb.*eber +dUdknkp.*eper;  dVdkrkn=dVdkn2 .*ener +dVdknkb.*eber +dVdknkp.*eper;
 dUdkrkb=dUdknkb.*ener +dUdkb2 .*eber +dUdkbkp.*eper;  dVdkrkb=dVdknkb.*ener +dVdkb2 .*eber +dVdkbkp.*eper;
 dUdkrkp=dUdknkp.*ener +dUdkbkp.*eber +dUdkp2 .*eper;  dVdkrkp=dVdknkp.*ener +dVdkbkp.*eber +dVdkp2 .*eper;
 %dUdkfkn=dUdkn2 .*enef +dUdknkb.*ebef +dUdknkp.*epef;  dVdkfkn=dVdkn2 .*enef +dVdknkb.*ebef +dVdknkp.*epef;
 %dUdkfkb=dUdknkb.*enef +dUdkb2 .*ebef +dUdkbkp.*epef;  dVdkfkb=dVdknkb.*enef +dVdkb2 .*ebef +dVdkbkp.*epef;
 %dUdkfkp=dUdknkp.*enef +dUdkbkp.*ebef +dUdkp2 .*epef;  dVdkfkp=dVdknkp.*enef +dVdkbkp.*ebef +dVdkp2 .*epef;
 dUdkzkn=dUdkn2 .*enez +dUdknkb.*ebez +dUdknkp.*epez;  dVdkzkn=dVdkn2 .*enez +dVdknkb.*ebez +dVdknkp.*epez;
 dUdkzkb=dUdknkb.*enez +dUdkb2 .*ebez +dUdkbkp.*epez;  dVdkzkb=dVdknkb.*enez +dVdkb2 .*ebez +dVdkbkp.*epez;
 dUdkzkp=dUdknkp.*enez +dUdkbkp.*ebez +dUdkp2 .*epez;  dVdkzkp=dVdknkp.*enez +dVdkbkp.*ebez +dVdkp2 .*epez;
 dUdr2= dUdrs .*dsdr +dUdrt .*dtdr;                    dVdr2= dVdrs .*dsdr +dVdrt .*dtdr;
 dUdz2= dUdzs .*dsdz +dUdzt .*dtdz;                    dVdz2= dVdzs .*dsdz +dVdzt .*dtdz;
 dUdrz= dUdrs .*dsdz +dUdrt .*dtdz;                    dVdrz= dVdrs .*dsdz +dVdrt .*dtdz;
 dUdrkr=dUdskr.*dsdr +dUdtkr.*dtdr;                    dVdrkr=dVdskr.*dsdr +dVdtkr.*dtdr;
 %dUdrkf=dUdskf.*dsdr +dUdtkf.*dtdr;                    dVdrkf=dVdskf.*dsdr +dVdtkf.*dtdr;
 dUdrkz=dUdskz.*dsdr +dUdtkz.*dtdr;                    dVdrkz=dVdskz.*dsdr +dVdtkz.*dtdr;
 dUdzkr=dUdskr.*dsdz +dUdtkr.*dtdz;                    dVdzkr=dVdskr.*dsdz +dVdtkr.*dtdz;
 %dUdzkf=dUdskf.*dsdz +dUdtkf.*dtdz;                    dVdzkf=dVdskf.*dsdz +dVdtkf.*dtdz;
 dUdzkz=dUdskz.*dsdz +dUdtkz.*dtdz;                    dVdzkz=dVdskz.*dsdz +dVdtkz.*dtdz;
 dUdkr2= dUdkrkn.*ener +dUdkrkb.*eber +dUdkrkp.*eper;  dVdkr2= dVdkrkn.*ener +dVdkrkb.*eber +dVdkrkp.*eper;
 %dUdkf2= dUdkfkn.*enef +dUdkfkb.*ebef +dUdkfkp.*epef;  dVdkf2= dVdkfkn.*enef +dVdkfkb.*ebef +dVdkfkp.*epef;
 dUdkz2= dUdkzkn.*enez +dUdkzkb.*ebez +dUdkzkp.*epez;  dVdkz2= dVdkzkn.*enez +dVdkzkb.*ebez +dVdkzkp.*epez;
 %dUdkrkf=dUdkrkn.*enef +dUdkrkb.*ebef +dUdkrkp.*epef;  dVdkrkf=dVdkrkn.*enef +dVdkrkb.*ebef +dVdkrkp.*epef;
 dUdkrkz=dUdkrkn.*enez +dUdkrkb.*ebez +dUdkrkp.*epez;  dVdkrkz=dVdkrkn.*enez +dVdkrkb.*ebez +dVdkrkp.*epez;
 %dUdkfkz=dUdkfkn.*enez +dUdkfkb.*ebez +dUdkfkp.*epez;  dVdkfkz=dVdkfkn.*enez +dVdkfkb.*ebez +dVdkfkp.*epez;
 dTdr2=  dVdr2  +pm.*dm1o2.*(-dm1o1.*(V.*dVdr -0.5*dUdr ).^2                   +dVdr.^2     +V.*dVdr2  -0.5*dUdr2  );
 dTdz2=  dVdz2  +pm.*dm1o2.*(-dm1o1.*(V.*dVdz -0.5*dUdz ).^2                   +dVdz.^2     +V.*dVdz2  -0.5*dUdz2  );
 dTdrz=  dVdrz  +pm.*dm1o2.*(-dm1o1.*(V.*dVdr -0.5*dUdr ).*(V.*dVdz -0.5*dUdz) +dVdr .*dVdz +V.*dVdrz  -0.5*dUdrz  );
 dTdrkr= dVdrkr +pm.*dm1o2.*(-dm1o1.*(V.*dVdr -0.5*dUdr ).*(V.*dVdkr-0.5*dUdkr)+dVdr .*dVdkr+V.*dVdrkr -0.5*dUdrkr );
 %dTdrkf= dVdrkf +pm.*dm1o2.*(-dm1o1.*(V.*dVdr -0.5*dUdr ).*(V.*dVdkf-0.5*dUdkf)+dVdr .*dVdkf+V.*dVdrkf -0.5*dUdrkf );
 dTdrkz= dVdrkz +pm.*dm1o2.*(-dm1o1.*(V.*dVdr -0.5*dUdr ).*(V.*dVdkz-0.5*dUdkz)+dVdr .*dVdkz+V.*dVdrkz -0.5*dUdrkz );
 dTdzkr= dVdzkr +pm.*dm1o2.*(-dm1o1.*(V.*dVdz -0.5*dUdz ).*(V.*dVdkr-0.5*dUdkr)+dVdz .*dVdkr+V.*dVdzkr -0.5*dUdzkr );
 %dTdzkf= dVdzkf +pm.*dm1o2.*(-dm1o1.*(V.*dVdz -0.5*dUdz ).*(V.*dVdkf-0.5*dUdkf)+dVdz .*dVdkf+V.*dVdzkf -0.5*dUdzkf );
 dTdzkz= dVdzkz +pm.*dm1o2.*(-dm1o1.*(V.*dVdz -0.5*dUdz ).*(V.*dVdkz-0.5*dUdkz)+dVdz .*dVdkz+V.*dVdzkz -0.5*dUdzkz );
 dTdkr2= dVdkr2 +pm.*dm1o2.*(-dm1o1.*(V.*dVdkr-0.5*dUdkr).^2                   +dVdkr.^2    +V.*dVdkr2 -0.5*dUdkr2 );
 %dTdkf2= dVdkf2 +pm.*dm1o2.*(-dm1o1.*(V.*dVdkf-0.5*dUdkf).^2                   +dVdkf.^2    +V.*dVdkf2 -0.5*dUdkf2 );
 dTdkz2= dVdkz2 +pm.*dm1o2.*(-dm1o1.*(V.*dVdkz-0.5*dUdkz).^2                   +dVdkz.^2    +V.*dVdkz2 -0.5*dUdkz2 );
 %dTdkrkf=dVdkrkf+pm.*dm1o2.*(-dm1o1.*(V.*dVdkr-0.5*dUdkr).*(V.*dVdkf-0.5*dUdkf)+dVdkr.*dVdkf+V.*dVdkrkf-0.5*dUdkrkf);
 dTdkrkz=dVdkrkz+pm.*dm1o2.*(-dm1o1.*(V.*dVdkr-0.5*dUdkr).*(V.*dVdkz-0.5*dUdkz)+dVdkr.*dVdkz+V.*dVdkrkz-0.5*dUdkrkz);
 %dTdkfkz=dVdkfkz+pm.*dm1o2.*(-dm1o1.*(V.*dVdkf-0.5*dUdkf).*(V.*dVdkz-0.5*dUdkz)+dVdkf.*dVdkz+V.*dVdkfkz-0.5*dUdkfkz);
end %if(Con,Amp)

if (rays.odeDim>=9 && evalAll)                     % --- Amplitude terms

%polcor pol=[]; dim1=size(DD,1); dim2=sqrt(size(DD,2));     % Polarization
%polcor for k=1:dim1
%polcor   dd=reshape(DD(k,:),dim2,dim2);
%polcor   new= ones(dim2,1); old=zeros(dim2,1); tol=1E-4; itmax=30; it=0;
%polcor   while (norm((new-old)./new,inf)>tol & it<=itmax)
%polcor     old=new; it=it+1; new=dd*new; new=new/norm(new);
%polcor   end
%polcor   if (it>itmax)
%polcor     disp(['dispertok: polarization not converged'])
%polcor   end
%polcor   pol=[pol new];
%polcor end;
                                            % Evolve in x-space
 dWr2Nx=dTdr2 +dTdrkr.*dWr2       +dTdrkz .*dWrz...
              +dTdrkr.*dWr2       +dTdrkz .*dWrz...
              +dTdkr2.*dWr2.*dWr2 +dTdkrkz.*dWr2.*dWrz...
              +dTdkz2.*dWrz.*dWrz +dTdkrkz.*dWrz.*dWr2;
 dWrzNx=dTdrz +dTdzkr.*dWr2       +dTdzkz .*dWrz...
              +dTdrkr.*dWrz       +dTdrkz .*dWz2...
              +dTdkr2.*dWr2.*dWrz +dTdkrkz.*dWr2.*dWz2...
              +dTdkz2.*dWrz.*dWz2 +dTdkrkz.*dWrz.*dWrz;
 dWz2Nx=dTdz2 +dTdzkr.*dWrz       +dTdzkz .*dWz2...
              +dTdzkr.*dWrz       +dTdzkz .*dWz2...
              +dTdkr2.*dWrz.*dWrz +dTdkrkz.*dWrz.*dWz2...
              +dTdkz2.*dWz2.*dWz2 +dTdkrkz.*dWz2.*dWrz;
 dlnE2Nx=(dTdrkr+dTdzkz +dTdkr2.*dWr2 +2.*dTdkrkz.*dWrz +dTdkz2.*dWz2);
 dphNx=-(kr.*dTdkr +kf.*dTdkf +kz.*dTdkz);
%polcor dphNx1=i*sum(conj(reshape(pol,dim1,dim2)) ...
%polcor                .* reshape(poldot,dim1,dim2),2)

                                            % Evolve in k-space
 dWr2Nk=-(dTdkr2 +dTdrkr.*dWr2       +dTdzkr.*dWrz...
                 +dTdrkr.*dWr2       +dTdzkr.*dWrz...
                 +dTdr2 .*dWr2.*dWr2 +dTdrz .*dWr2.*dWrz...
                 +dTdz2 .*dWrz.*dWrz +dTdrz .*dWrz.*dWr2);
 dWrzNk=-(dTdkrkz+dTdrkz.*dWr2       +dTdzkz.*dWrz...
                 +dTdrkr.*dWrz       +dTdzkr.*dWz2...
                 +dTdr2 .*dWr2.*dWrz +dTdrz .*dWr2.*dWz2...
                 +dTdz2 .*dWrz.*dWz2 +dTdrz .*dWrz.*dWrz);
 dWz2Nk=-(dTdkz2 +dTdrkz.*dWrz       +dTdzkz.*dWz2...
                 +dTdrkz.*dWrz       +dTdzkz.*dWz2...
                 +dTdr2 .*dWrz.*dWrz +dTdrz .*dWrz.*dWz2...
                 +dTdz2 .*dWz2.*dWz2 +dTdrz .*dWz2.*dWrz);
 dlnE2Nk=-(dTdrkr+dTdzkz +dTdr2.*dWr2 +2.*dTdrz.*dWrz +dTdz2.*dWz2);
 dphNk=-(rr.*dTdr +zz.*dTdz);
%polcor dphNx1=i*sum(conj(reshape(pol,dim1,dim2)) ...
%polcor                .* reshape(poldot,dim1,dim2),2);

 % Choose which space [I think this needs to be modified to use only the
 % index of inKspace which corresponds the the ray we are calculating.]
 dWr2N =rays.inKspace'.*dWr2Nk +(1.-rays.inKspace').*dWr2Nx;
 dWrzN =rays.inKspace'.*dWrzNk +(1.-rays.inKspace').*dWrzNx;
 dWz2N =rays.inKspace'.*dWz2Nk +(1.-rays.inKspace').*dWz2Nx;
 dphN  =rays.inKspace'.*dphNk  +(1.-rays.inKspace').*dphNx;
 dlnE2N=rays.inKspace'.*dlnE2Nk+(1.-rays.inKspace').*dlnE2Nx;

 Qdep=zeros(size(plasma.acharge));                 % Power dissipation
 one=ones(size(plasma.amass'))';
 lrmr2=vth2./omc2; kvt=abs(kp)*one.*vth;
 En=pol3(:,1); argp0=(plasma.omega      )./kvt;
 Eb=pol3(:,2); argp1=(plasma.omega+1*omc)./kvt; argm1=(plasma.omega-1*omc)./kvt;
 Ep=pol3(:,3); argp2=(plasma.omega+2*omc)./kvt; argm2=(plasma.omega-2*omc)./kvt;
 f1=2*plasma.omega*omc./(kp*one.*vth2);
 iomBp=i*(kn.*Eb-kb.*En)*one;   TTMP2=conj(iomBp).*iomBp;
 land =f1.*(Ep*one)-iomBp;      land2=conj(land) .*land;
 kpls=i*kn-kb; kmns=i*kn+kb;
 Em1=En-i*Eb;                   Em1sq=conj(Em1)  .*Em1*one;
 Ep1=En+i*Eb;                   Ep1sq=conj(Ep1)  .*Ep1*one;
 Em2=kmns.*Em1;                 Em2sq=conj(Em2)  .*Em2*one;
 Ep2=kpls.*Ep1;                 Ep2sq=conj(Ep2)  .*Ep2*one;
                                             % Interactions with species
 Qresn=   exp(-argm1.^2).* Em1sq + ...       %   1st cyclotron harmonic
   lrmr2.*exp(-argm2.^2).* Em2sq + ...       %   2nd cyclotron harmonic
   lrmr2.*exp(-argp0.^2).*(TTMP2+land2 );    %   Landau damping & magn pumping
 Qres =sqrt(pi)./(plasma.omega*abs(kp)*one).* ...   % Damping resonant interactions
       omp2./vth.*Qresn;                     %   (8*pi/omega)*P -> JanRev 6.21
 Qvisc=(vth2.*(N2*one/coomsq)./omc2).^2;     % Anomalous viscosity
 Qvis =0E-3*(plasma.omega./nu).*Qvisc;
 Qtot =Qres+Qvis;                            % Total damping
 dlnE2N=dlnE2N +2*sum(Qtot,2);               % !sign(dTdom)=-1! -> JanRev 7.36
 Pdep=-Qtot.*(exp(lnE2)*one)*2*plasma.omega;        % Power deposition per species
%
%%%%%%%%%%
% Strange: Em1->0 at cyclotron resonance and not Ep1.
%          Negative magn field or sign in time evolution?
%
if (sys.dbflag)                                  % Print damping rate
  damp=[Qres Qvisc]                          %   gamma/omega
end
%%%%%%%%%%
end; %Amp
%
% ===== Caustics ==============================================================
%
mon1=zeros(dimy,1); tune=4;
if (rays.odeDim>=9 && strcmp(oper,'Mon'))
 for k=1:dimy
   %mon1(k)=1.-tune*abs(ydot(k,8))/plasma.omega;    % derivative of log(E^2)
   mon1(k) = (~rays.inKspace(1)*3000 + rays.inKspace(1)*0.01)...
       - norm(y(k,5:7));   % size of elements of focusing tensor
 end
end
%
% ===== Mode conversion ======================================================
%
% ----- Estimates based on expansion at given point --------------------------
if (strcmp(oper,'Mon')||strcmp(oper,'Mch')||strcmp(oper,'Sdl')||...
    strcmp(oper,'Cnv')||strcmp(oper,'Trs'))
  detH=H11.*H22-H12.^2;                     % Det and inverse of Hessian
  Hm11=H22./detH;  Hm12=-H12./detH;
  Hm21=Hm12;       Hm22= H11./detH;
  eta2=0.5./sqrt(abs(detH)) ...             % Coupling coef (Tracy, eq.26) 
     .*abs(Hm11.* dUdq.^2 +2.*Hm12.*(dUdq.*dUdp) +Hm22.*dUdp.^2);
  qst=(Hm11.*dUdq +Hm12.*dUdp)*ones(1,4);   % Converted ray (Tracy, eq.24)
  pst=(Hm21.*dUdq +Hm22.*dUdp)*ones(1,4);
  zinzst=-(qst.*eq+pst.*ep);                % Vector from z* to z0, eq.8,7)
  yg=yin+2*[zinzst zeros(dimy,rays.odeDim-4)];   % Guess transmitted ray
end;
%
% ----- Store coordinates of incoming ray ------------------------------------
if (strcmp(oper,'Mon'))
  rays.yalf0=yin;                                % Store for subsequent use
  rays.eta2est=eta2;
%
% ----- Newton iteration to compute saddle point
elseif (strcmp(oper,'Sdl'))
  y=yin+[zinzst zeros(dimy,rays.odeDim-4)];      % The saddle point z*
%
% ----- Matching to local form to compute Eikonal phase ----------------------
elseif (strcmp(oper,'Mch')||strcmp(oper,'Trs')||strcmp(oper,'Cnv'))
  J2=[0 1; -1 0];                           % Symplectic operators
  J4=[0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0];
  T0=T;                                     % Expansion about z*
  T1=[dUdr dUdz dUdkr dUdkz];               %   gradient in (x,k)
  T2=[dUdr2  dUdrz  dUdrkr  dUdrkz          %   hessian in (x,k)
      dUdrz  dUdz2  dUdzkr  dUdzkz
      dUdrkr dUdzkr dUdkr2  dUdkrkz
      dUdrkz dUdzkz dUdkrkz dUdkz2];
  [M,val]=eig(J4*T2);                       % Transformation to normal form
  [val,ind]=sort(real(diag(val)));          %   increasing order real part
  opposite=val(4)./val(1);
  separate=norm(val(4)/val(3));
  if (abs(1+opposite)<0.1 && separate>4)     % Well behaved hyperbola
    vp=M(:,ind(4)); vm=M(:,ind(1));         %   direction of asymptotes
    dim2=sqrt(size(DD,2)); cgD12=conj(gD12); 
    if ((strcmp(plasma.MODEL(4:6),'2x2')))
      gdv=[ gD11*vp gD12*vp cgD12*vp gD22*vp;
            gD11*vm gD12*vm cgD12*vm gD22*vm];
    elseif ((strcmp(plasma.MODEL(4:6),'3x3')))
      cgD13=conj(gD13); cgD23=conj(gD23); 
      gdv=[ gD11*vp  gD12*vp gD13*vp ... 
           cgD12*vp  gD22*vp gD23*vp ...
           cgD13*vp cgD23*vp gD23*vp;
            gD11*vm  gD12*vm gD13*vm ... 
           cgD12*vm  gD22*vm gD23*vm ...
           cgD13*vm cgD23*vm gD23*vm];
    end
    pol=[];
    for k=1:2                               %   uncoupled polarizations
      gd=reshape(gdv(k,:),dim2,dim2);
      new= ones(dim2,1); old=zeros(dim2,1); tol=1E-4; itmax=100; it=0;
      while (norm((new-old)./new,inf)>tol && it<=itmax)
        old=new; it=it+1; new=gd*new; new=new/norm(new);
      end
      pol=[pol new];
    end
    if ((strcmp(plasma.MODEL(4:6),'2x2')))         % Use trick to expand e*.D.e
      gD = [gD11' gD12' cgD12' gD22'];      %   create matrix and reshape
    elseif ((strcmp(plasma.MODEL(4:6),'3x3')))
      gD = [gD11' gD12' gD13' cgD12' gD22' gD23' cgD13' cgD23' gD33'];
    end
    gdalf = gD*reshape(conj(pol(:,1)*pol(:,1)'),1,size(gD,2))';
    gdlam = gD*reshape(conj(pol(:,2)*pol(:,2)'),1,size(gD,2))';
    braket=gdalf'*J4*gdlam;                 % Coupling coefficient
    eta=pol(:,1)'*reshape(DD,dim2,dim2)*pol(:,2)/sqrt(braket);
    eta2=eta*conj(eta);
    tau=exp(-pi*eta2);
    beta=sqrt(2*pi*tau)/(eta*cgamma(-i*eta2));

    if(rays.odeDim>=9)                           % Matching to local 2x2 form
     za=gdalf(1:2); zb=gdalf(3:4);
     zc=gdlam(1:2); zd=gdlam(3:4);
     Salf =[rays.yalf0(5:6); rays.yalf0(6:7)];        % Incoming data inconsistent!
                                            %   should use anterior time
     Slam0=(zc'*zd)/(zd'*zd)^2*zd*zd' ...   % Outgoing Hessian
          -(zc*zd' +zd*zc')/(zd'*zd);
     Slam1=(J2*zd)*(J2*zd)';
     zvc=J2*(zc+Salf*zd);
     lambda=-(za+Slam0*zb)'*zvc/(zb'*Slam1*zvc);
     Slam = Slam0+lambda*Slam1;
    end % Amp
  else
    tau=0
  end % hyperbola
%
% ----- Forget stored saddle point -------------------------------------------
else
  ysdl=[];
end
%
% ----- Lauch transmitted ray ------------------------------------------------
if (strcmp(oper,'Trs'))
  zinzst=[y(1:4)-rays.yalf0(1:4) ...   % Initialize transmitted in (x,k)
          zeros(dimy,rays.odeDim-4)];       %   direction z0 -> z*
  yg=rays.yalf0-2.2*zinzst;                 %   -2 is a guess (eq.28)
  fact=fzero('dispertok',2., ...            %   fact is the actual value
             optimset('Display','off'), ...
             rays.yalf0',zinzst',zeros(size(zinzst))',zeros(1,3),'Dsp');
  y=rays.yalf0+fact*zinzst;            %   exactly on dispersion manifold
  if (rays.odeDim>=9)                       % S-matrix remains unchanged
   y(8)=rays.y(8)-2*pi*eta2;           % Transmitted log(E^2) and phase
   %ampi=sqrt(exp(rays.yalf0(8)));           %   incoming data approximative,
   %phsi=rays.yalf0(9);                      %   take far from coupling region!
   y(10:end)=zeros(size(plasma.amass))';% Reset power deposition
%%%%%%
%  zbsq=zb'*zb;  zbp=J2*zb;                 % Only for checking purposes
%  zdsq=zd'*zd;  zdp=J2*zd;                 %  uncoupled phase function
%  zxi=yalf0(1:2)';                         %   before conversion
%  alpha=1./(zxi'*zbp)^2 * ...
%      (2*phsi -1/zbsq^2*(za'*zb )*(zb'*zxi)^2 ...
%              +2/zbsq  *(za'*zxi)*(zb'*zxi) );
%  zxo=y(1:2)';                             %   after conversion
%  phsalfo=0.5*(za'*zb )*(zb'*zxo)^2/zbsq^2 -(za'*zxo)*(zb'*zxo)/zbsq + ...
%          0.5*alpha*(zxo'*zbp);
%  phslamo=0.5*(zc'*zd )*(zd'*zxo)^2/zdsq^2 -(zc'*zxo)*(zd'*zxo)/zdsq + ...
%          0.5*lambda*(zxo'*zdp);
%                                           %  conversion line function H~
%  zx=zxi; gphs=(za'*zb)*(zb'*zx)*zb/zbsq^4 ...
%            -((zb'*zx)*za +(za'*zx)*zb)/zbsq +alpha *(zbp'*zx)*zbp;
%          htalfi=(za'*zx)+(zb'*gphs);
%          gphs=(zc'*zd)*(zd'*zx)*zd/zdsq^4 ...
%             -((zd'*zx)*zc +(zc'*zx)*zd)/zdsq +lambda*(zdp'*zx)*zdp;
%          htlami=(zc'*zx)+(zd'*gphs);
%  zx=zxo; gphs=(za'*zb)*(zb'*zx)*zb/zbsq^4 ...
%            -((zb'*zx)*za +(za'*zx)*zb)/zbsq +alpha *(zbp'*zx)*zbp;
%          htalfo=(za'*zx)+(zb'*gphs);
%          gphs=(zc'*zd)*(zd'*zx)*zd/zdsq^4 ...
%             -((zd'*zx)*zc +(zc'*zx)*zd)/zdsq +lambda*(zdp'*zx)*zdp;
%          htlamo=(zc'*zx)+(zd'*gphs);
%                                           %  solution 2x2 form
%  Ealfo=tau*exp(i*(phsalfo-phsi))*ampi * ...
%      (htlamo/htlami)^(i*conj(eta)*eta/braket);
%%%%%
  end; % Amp
%
% ----- Modify outgoing converted ray ----------------------------------------
elseif (strcmp(oper,'Cnv'))
  disp( sprintf('   |eta2|=%g   |est2|=%g',norm(eta2),rays.eta2est) )
  disp( sprintf('   tau=%g  |beta|=%g  angle=%g',tau,norm(beta),angle(beta)))
%
  % Save conversion data
  rays.tau = tau;   rays.beta = beta;
  
  rays.y=rays.yalf0;                        % Initialize z2out in (x,k)
  if (rays.odeDim>=9)
    disp( sprintf('   Salf = [%0.2g %0.2g; %0.2g %0.2g]',reshape(Salf,1,4)) )
    disp( sprintf('   Slam = [%0.2g %0.2g; %0.2g %0.2g]',reshape(Slam,1,4)) )
    nrmbeta=norm(beta);
    phsbeta=angle(beta);
    rays.y(5)=Slam(1,1);                         % Converted S-matrix
    rays.y(6)=0.5*(Slam(1,2)+Slam(2,1));
    rays.y(7)=Slam(2,2);
    rays.y(8)=rays.y(8)+log(nrmbeta^2);          % Converted log(E^2)
    rays.y(9)=rays.y(9)+phsbeta;                 %           Eikonal phase
  end
end
%
% ===== Variable output dimensions ===========================================
%
switch oper
case {'Frq'}                            % --- Analytical frequencies
  if (nspec==3)
    omii= sqrt(omp2(:,2).*(omc(:,3).^2)+omp2(:,3).*(omc(:,2).^2))./ ...
          sqrt(omp2(:,2)+omp2(:,3));
    param=cat(2,omii,omc(:,2),omc(:,3))./(2*pi);
  else
    param=omc(:,2:nspec)./(2*pi);
  end;
  out = cat(1, param);
case {'Dsp','Msw'}                      % --- Dispersion relations
%  out = real(prod(omc2Mom2,2).*U);
  out = real(U);
case {'Ten'}                            % --- Dispersion tensor
  out = DD;
case {'Pol'}                            % --- Polarization of rays
  out = pol;
case {'Ant'}                            % --- IC for a wavefront aligned
  dzdr=-dsdr./dsdz;                     %     flux surfaces
  % NB: the above is bad if starting exactly on the midplane
  %     since there dsdz=0;
  matrx=[ [1 2*dzdr dzdr.^2]; ...
          [dTdkr dTdkz 0]; ...
          [0 dTdkr dTdkz]   ];
  rhs=[0; -dTdr; -dTdz];  
  dW=matrx\rhs; 
  dWr2=dW(1); dWrz=dW(2); dWz2=dW(3); 
  lnE2=zeros(size(dWr2));phs=zeros(size(dWr2));
  out=real([dWr2,dWrz,dWz2,lnE2,phs]');
case {'Trj'}                            % --- ODE trajectory
%  rays.dUdomHist = [rays.dUdomHist dUdom];
%  dUdomM=ones(rays.odeDim,1)*(1./dUdom)';  % fn of t, not sigma
  dUdomM=ones(rays.odeDim,1);           % solve as fn of sigma, not t
  out=real(reshape(dUdomM.*cat(1,-dUdkr',-dUdkz',dUdr',dUdz'),dimytot,1));

  % Use eigenvalue as ray hamiltonian instead of det.
  % This is wrong if using the 3x3 model, use dUd... above
  %dTdomM=ones(rays.odeDim,1)*(1./dTdom)'; % -5.8611e-007
%   dTdomM=-ones(rays.odeDim,1); % use ray parameter, not time
%   out=real(reshape(dTdomM.*cat(1,-dTdkr',-dTdkz',dTdr',dTdz'),dimytot,1));

case {'Eig'}                            % --- Eigenvalue
  out=eig2;
case {'Mon'}                            % --- Conversion monitors
  tmp=real(cat(1,mon1',mon2',eta2',yg(:,1)',yg(:,3)'));
  out=reshape(tmp,1,numel(tmp));
case {'Mch'}                            % --- Conversion monitors
  tmp=real(cat(1,tau'));
  out=reshape(tmp,1,numel(tmp));
case {'Sdl'}                            % --- Conversion find saddle point
  out=[y(1:4) zeros(1,rays.odeDim-4)];
case {'Cnv','Trs'}                      % --- Quadratic conversion 
  tmp=cat(2,y);
  out=real(reshape(tmp',numel(tmp),1));
case {'Amp'}                            % --- ODE trajectory & amplitude
%  dTdomM=ones(rays.odeDim,1)*(1./dTdom)';    % !!Note that this is negative!!
  dTdomM=-ones(rays.odeDim,1); % use ray parameter, not time
  out=reshape(dTdomM.*cat(1, ...        % ==> Invert all the signs in RHS
      -dTdkr',-dTdkz',dTdr',dTdz', ...
      dWr2N(1)',dWrzN(1)',dWz2N(1)',dlnE2N(1)',dphN(1)', ...
      Pdep'),dimytot,1);
%      dWr2N',dWrzN',dWz2N',dlnE2N',dphN', ...  these are too big if
%      inKspace > 1, so use dWr2N(1)', etc. (mod by Steve)
  out=real(out);
case {'Sgn'}
  %out = sign(dUdom); % return the sign of the derivative of time wrt ray parameter
  out = (dUdom); % return the sign of the derivative of time wrt ray parameter
otherwise
  error(['Unknown oper ''' oper '''.']);
end

function out= disp_eig(yv,oper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DISPERSION --  Plasma dispersion characteristics
%                 Part of the RAYcON package
%
%  yin     Position in phase space [r,z,k_r,k_z]
%  oper    
%          'Msw' Dispersion relation approximated 1x1
%          'Dsp' Disperison relation according to MODEL
%          'Ant' Antenna IC
%          'Trj' RHS for trajectory
%          'Sgn' Sign of dt/dsigma, set by dUdom
%          'Pol' Polarization vector, eigenvector of D matrix
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
global cnst 
global plasma rays

i = complex(0,1);
isCmplx=0;
%
% ----- Vectorized input quantities ------------------------------------------
%
if numel(yv) == 4
    yin = yv.';
    sym_mat = eye(4);
    eval2nd = 0;
else
    yin = yv;
    % now assume one ray, and tracing position and symplectic matrix
    sym_mat = reshape(yin(5:20),4,4);
    yin = yin(1:4).';
    eval2nd = 1;
end

if strcmp(oper,'Trj')
    eval1st = 1;
else
    eval1st = 0;
    eval2nd = 0;
end

%
% ----- Adjusted position & wave vector in cylindrical coordinates -----------
zro=zeros(size(yin(:,1))); one=ones(size(zro));
if (rays.odeDim==4||rays.odeDim>=9)
 rr=yin(:,1); zz=yin(:,2); kr=yin(:,3); kz=yin(:,4);
 kf=one.*plasma.kant(2);                            % fi=omega*tim/(kf*rr);
end
rho=sqrt((rr-plasma.r0).^2+zz.^2);
theta=atan2(plasma.r0-rr,-zz)+pi/2;
%
%
% ===== Local plasma parameters ==============================================
%
% ----- Magnetic field and topology ------------------------------------------
  [b,dbds,dbdt,bp,sflx,dsdr,dsdz,dtdr,dtdz, ...
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
         depezdr2, depezdrz, depezdz2]=magnetic(rho,theta);



%
% ----- Wave vector and refraction index -------------------------------------
coom=cnst.c/plasma.omega; coomsq=coom^2; om2=plasma.omega^2;
kn=kr.*ener +kf.*enef +kz.*enez;  Nn=coom*kn;  Nn2=Nn.^2;
kb=kr.*eber +kf.*ebef +kz.*ebez;  Nb=coom*kb;  Nb2=Nb.^2;
kp=kr.*eper +kf.*epef +kz.*epez;  Np=coom*kp;  Np2=Np.^2; %N2=Nn2+Nb2+Np2;

Nr = coom*kr;
Nf = coom*kf;
Nz = coom*kz;

%
% ----- Parabolic profiles and logarithmic derivatives -----------------------
nspec=size(plasma.amass,2);
p=(1-sflx.^2*plasma.na);
dLNpds = zeros(length(sflx),nspec);
dLNpds2 = dLNpds;
n = dLNpds;
T = dLNpds;
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
% second derivatives
dLNomcds2=(dbds2./b-(dbds.^2   )./(b.^2))*one;
dLNomcdst=(dbdst./b-(dbds.*dbdt)./(b.^2))*one;
dLNomcdt2=(dbdt2./b-(dbdt.^2   )./(b.^2))*one;

switch plasma.MODEL(1:6)
%
% +++++ Cold plasma model ++++++++++++++++++++++++++++++++++++++++++++++++++++
%
case {'cld2x2'}

  % --- Elementary functions and derivatives ---------------------------------
  if (isCmplx && ~eval1st)
    iomceps=0.003*i*om2;
  else
    iomceps=0;
  end
  omc2Mom2= omc2-om2-2*iomceps;              % S,D,P
  Si=omp2./omc2Mom2;   S=1+sum( Si ,2);
  Di=(omc/plasma.omega).*Si;  D=  sum( Di ,2);
  % Pi=omp2/om2;         P=1-sum( Pi ,2);

                                              % 1st derivatives S,D,P
  zro=zeros(size(S));
  dLNSids  =dLNnds-2.*omc2./omc2Mom2.*dLNomcds; dSds = sum( Si.*dLNSids ,2);
  dLNSidt  =      -2.*omc2./omc2Mom2.*dLNomcdt; dSdt = sum( Si.*dLNSidt ,2);
  dLNSidom =2*plasma.omega./omc2Mom2;           dSdom= sum( Si.*dLNSidom,2);
  dLNDids  =dLNSids + dLNomcds;                 dDds = sum( Di.*dLNDids ,2);
  dLNDidt  =dLNSidt + dLNomcdt;                 dDdt = sum( Di.*dLNDidt ,2);
  dLNDidom =(3*om2-omc2)./(plasma.omega*omc2Mom2);     dDdom= sum( Di.*dLNDidom,2);

  % % These are only needed in the 3x3 model
  % zroi=zeros(size(Si));                         
  % dLNPids =dLNnds;                             
  % dPds =-sum( Pi.*dLNPids ,2);
  % dLNPidt =zroi;                               
  % dPdt = zro;
  % dLNPidom=2/plasma.omega;                     
  % dPdom= sum( Pi.*dLNPidom,2);
  
  % second derivatives
  if eval2nd
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

%    % These are only needed in the 3x3 model
%    dLNPids2= dLNnds2;   dLNPidst= zroi;  dLNPidt2= zroi;
%    dPids2= Pi.*(dLNPids2 +dLNPids.*dLNPids);   dPds2=sum( dPids2 ,2);
%    dPidst= Pi.*(dLNPidst +dLNPids.*dLNPidt);   dPdst=sum( dPidst ,2);
%    dPidt2= Pi.*(dLNPidt2 +dLNPidt.*dLNPidt);   dPdt2=sum( dPidt2 ,2);
  end
  
  % --- Tensor elements and derivatives --------------------------------------

    D11 = Nb2+Np2-S;
    D12 =-Nn.*Nb-D*i;
    D22 = Nn2+Np2-S;
    
    cD12 = conj(D12);
    
    % Dispersion tensor
    DD=[D11 D12 cD12 D22];
    % We will assume that the input yv is size [1,4], then reshape works
    DD=reshape(DD,2,2);
    [evects,evals]=eig(DD);           % Two values
    evs = diag(evals);
    ind = find(abs(evs)==min(abs(evs)));
    U = evs(ind); % The eval nearest zero is the one we want
    pol = evects(:,ind); % get the column with the right evector
    
    % Conversion monitor
    mon2 = abs(D11 + D22); % trace DD

    if eval1st
        dD11ds =-dSds;          dD11dt =-dSdt;
        dD12ds =-dDds*i;        dD12dt =-dDdt*i;
        dD22ds =-dSds;          dD22dt =-dSdt;

        dD11dkn= zro;           dD11dkb= Nb*coom*2;     dD11dkp= Np*coom*2;
        dD12dkn=-Nb*coom;       dD12dkb=-Nn*coom;       dD12dkp= zro;
        dD22dkn= Nn*coom*2;     dD22dkb= zro;           dD22dkp= Np*coom*2;

        dD11dkr=dD11dkn.*ener +dD11dkb.*eber +dD11dkp.*eper;
        dD11dkz=dD11dkn.*enez +dD11dkb.*ebez +dD11dkp.*epez;
        dD12dkr=dD12dkn.*ener +dD12dkb.*eber +dD12dkp.*eper;
        dD12dkz=dD12dkn.*enez +dD12dkb.*ebez +dD12dkp.*epez;
        dD22dkr=dD22dkn.*ener +dD22dkb.*eber +dD22dkp.*eper;
        dD22dkz=dD22dkn.*enez +dD22dkb.*ebez +dD22dkp.*epez;

        dD11dom=-2/plasma.omega*(Nb2+Np2)-dSdom;
        dD12dom=-2/plasma.omega*(Nn.*Nb) -dDdom*i;
        dD22dom=-2/plasma.omega*(Nn2+Np2)-dSdom;

        % Sub-blocks of (dispersion matrix - eigenval*id)
        subsum = D11+D22-2*U;

        % Convert to r, z coordinates.
        %  Note that there are extra correction terms due to field
        %  curvature.  They are derivatives of the trasformation elements:
        %  (d/dr)(ebef) etc.
        
        dNndr = Nr.*denerdr + Nf.*denefdr + Nz.*denezdr;
        dNndz = Nr.*denerdz + Nf.*denefdz + Nz.*denezdz;
        dNbdr = Nr.*deberdr + Nf.*debefdr + Nz.*debezdr;
        dNbdz = Nr.*deberdz + Nf.*debefdz + Nz.*debezdz;
        dNpdr = Nr.*deperdr + Nf.*depefdr + Nz.*depezdr;
        dNpdz = Nr.*deperdz + Nf.*depefdz + Nz.*depezdz;
        
        dD11dr = dD11ds.*dsdr + dD11dt.*dtdr + 2*(Nb.*dNbdr + Np.*dNpdr);
        dD11dz = dD11ds.*dsdz + dD11dt.*dtdz + 2*(Nb.*dNbdz + Np.*dNpdz);
        dD12dr = dD12ds.*dsdr + dD12dt.*dtdr - Nn.*dNbdr - Nb.*dNndr;
        dD12dz = dD12ds.*dsdz + dD12dt.*dtdz - Nn.*dNbdz - Nb.*dNndz;
        dD22dr = dD22ds.*dsdr + dD22dt.*dtdr + 2*(Nn.*dNndr + Np.*dNpdr);
        dD22dz = dD22ds.*dsdz + dD22dt.*dtdz + 2*(Nn.*dNndz + Np.*dNpdz);

        % Gather things together to do all derivatives at the same time
        dD11_vec = [dD11dom dD11dr dD11dz dD11dkr dD11dkz];
        dD12_vec = [dD12dom dD12dr dD12dz dD12dkr dD12dkz];
        dD22_vec = [dD22dom dD22dr dD22dz dD22dkr dD22dkz];

        % Derivatives of the eigenvalue
        %dU_vec   = [dUdom dUds dUdt dUdkn dUdkb dUdkp];
        dU_vec = (dD11_vec.*(D22-U) +dD22_vec.*(D11-U))./subsum...
               - 2*real(cD12.*dD12_vec)./subsum;

    end
    
    % terms needed for second derivative of hamiltonian
    if eval2nd
        
        dSdr2 = dsdr2.*dSds + dtdr2.*dSdt ...
            + dsdr.*(dsdr.*dSds2 + dtdr.*dSdst) + dtdr.*(dsdr.*dSdst + dtdr.*dSdt2);
        dSdrz = dsdrz.*dSds + dtdrz.*dSdt ...
            + dsdz.*(dsdr.*dSds2 + dtdr.*dSdst) + dtdz.*(dsdr.*dSdst + dtdr.*dSdt2);
        dSdz2 = dsdz2.*dSds + dtdz2.*dSdt ...
            + dsdz.*(dsdz.*dSds2 + dtdz.*dSdst) + dtdz.*(dsdz.*dSdst + dtdz.*dSdt2);
        dDdr2 = dsdr2.*dDds + dtdr2.*dDdt ...
            + dsdr.*(dsdr.*dDds2 + dtdr.*dDdst) + dtdr.*(dsdr.*dDdst + dtdr.*dDdt2);
        dDdrz = dsdrz.*dDds + dtdrz.*dDdt ...
            + dsdz.*(dsdr.*dDds2 + dtdr.*dDdst) + dtdz.*(dsdr.*dDdst + dtdr.*dDdt2);
        dDdz2 = dsdz2.*dDds + dtdz2.*dDdt ...
            + dsdz.*(dsdz.*dDds2 + dtdz.*dDdst) + dtdz.*(dsdz.*dDdst + dtdz.*dDdt2);

        dNndr2 = Nr.*denerdr2 + Nf.*denefdr2 + Nz.*denezdr2;
        dNndrz = Nr.*denerdrz + Nf.*denefdrz + Nz.*denezdrz;
        dNndz2 = Nr.*denerdz2 + Nf.*denefdz2 + Nz.*denezdz2;
        dNbdr2 = Nr.*deberdr2 + Nf.*debefdr2 + Nz.*debezdr2;
        dNbdrz = Nr.*deberdrz + Nf.*debefdrz + Nz.*debezdrz;
        dNbdz2 = Nr.*deberdz2 + Nf.*debefdz2 + Nz.*debezdz2;
        dNpdr2 = Nr.*deperdr2 + Nf.*depefdr2 + Nz.*depezdr2;
        dNpdrz = Nr.*deperdrz + Nf.*depefdrz + Nz.*depezdrz;
        dNpdz2 = Nr.*deperdz2 + Nf.*depefdz2 + Nz.*depezdz2;
        
        dD11dr2=2*(dNbdr.*dNbdr + Nb.*dNbdr2 + dNpdr.*dNpdr + Np.*dNpdr2) - dSdr2;
        dD11drz=2*(dNbdr.*dNbdz + Nb.*dNbdrz + dNpdr.*dNpdz + Np.*dNpdrz) - dSdrz;        
        dD11dz2=2*(dNbdz.*dNbdz + Nb.*dNbdz2 + dNpdz.*dNpdz + Np.*dNpdz2) - dSdz2;
        dD12dr2=-(Nn.*dNbdr2 + 2*dNndr.*dNbdr + Nb.*dNndr2)-i*dDdr2;
        dD12drz=-(Nn.*dNbdrz + dNndr.*dNbdz + dNndz.*dNbdr + Nb.*dNndrz)-i*dDdrz;
        dD12dz2=-(Nn.*dNbdz2 + 2*dNndz.*dNbdz + Nb.*dNndz2)-i*dDdz2;
        dD22dr2=2*(dNndr.*dNndr + Nn.*dNndr2 + dNpdr.*dNpdr + Np.*dNpdr2) - dSdr2;
        dD22drz=2*(dNndr.*dNndz + Nn.*dNndrz + dNpdr.*dNpdz + Np.*dNpdrz) - dSdrz;
        dD22dz2=2*(dNndz.*dNndz + Nn.*dNndz2 + dNpdz.*dNpdz + Np.*dNpdz2) - dSdz2;
        
        dD11drkr = 2*coom*(eber.*dNbdr + Nb.*deberdr + eper.*dNpdr + Np.*deperdr);
        dD11drkz = 2*coom*(ebez.*dNbdr + Nb.*debezdr + epez.*dNpdr + Np.*depezdr);
        dD11dzkr = 2*coom*(eber.*dNbdz + Nb.*deberdz + eper.*dNpdz + Np.*deperdz);
        dD11dzkz = 2*coom*(ebez.*dNbdz + Nb.*debezdz + epez.*dNpdz + Np.*depezdz);
        
        dD12drkr = -coom*(denerdr.*Nb + dNndr.*eber + ener.*dNbdr + Nn.*deberdr);
        dD12drkz = -coom*(denezdr.*Nb + dNndr.*ebez + enez.*dNbdr + Nn.*debezdr);
        dD12dzkr = -coom*(denerdz.*Nb + dNndz.*eber + ener.*dNbdz + Nn.*deberdz);
        dD12dzkz = -coom*(denezdz.*Nb + dNndz.*ebez + enez.*dNbdz + Nn.*debezdz);
        
        dD22drkr = 2*coom*(ener.*dNndr + Nn.*denerdr + eper.*dNpdr + Np.*deperdr);
        dD22drkz = 2*coom*(enez.*dNndr + Nn.*denezdr + epez.*dNpdr + Np.*depezdr);
        dD22dzkr = 2*coom*(ener.*dNndz + Nn.*denerdz + eper.*dNpdz + Np.*deperdz);
        dD22dzkz = 2*coom*(enez.*dNndz + Nn.*denezdz + epez.*dNpdz + Np.*depezdz);


        dD11dkr2 = 2*coomsq*(eber.^2 + eper.^2);
        dD11dkz2 = 2*coomsq*(ebez.^2 + epez.^2);
        dD11dkrkz= 2*coomsq*(eber.*ebez + eper.*epez);
        dD12dkr2 = -coomsq*(ener.*eber + ener.*eber);
        dD12dkz2 = -coomsq*(ener.*ebez + enez.*eber);
        dD12dkrkz= -coomsq*(enez.*ebez + enez.*ebez);
        dD22dkr2 = 2*coomsq*(ener.^2 + eper.^2);
        dD22dkz2 = 2*coomsq*(enez.^2 + epez.^2);
        dD22dkrkz= 2*coomsq*(ener.*enez + eper.*epez);

        % dD11_mat = [dr dz dkr dkz].'*[dr dz dkr dkz];
        dD11_mat = [dD11dr2   dD11drz   dD11drkr  dD11drkz;
                    dD11drz   dD11dz2   dD11dzkr  dD11dzkz;
                    dD11drkr  dD11dzkr  dD11dkr2  dD11dkrkz;
                    dD11drkz  dD11dzkz  dD11dkrkz dD11dkz2 ];
        dD12_mat = [dD12dr2   dD12drz   dD12drkr  dD12drkz;
                    dD12drz   dD12dz2   dD12dzkr  dD12dzkz;
                    dD12drkr  dD12dzkr  dD12dkr2  dD12dkrkz;
                    dD12drkz  dD12dzkz  dD12dkrkz dD12dkz2 ];
        dD22_mat = [dD22dr2   dD22drz   dD22drkr  dD22drkz;
                    dD22drz   dD22dz2   dD22dzkr  dD22dzkz;
                    dD22drkr  dD22dzkr  dD22dkr2  dD22dkrkz;
                    dD22drkz  dD22dzkz  dD22dkrkz dD22dkz2 ];

        % for the second order derivatives, we don't need dom
        dD11_vec = [dD11dr dD11dz dD11dkr dD11dkz];
        dD12_vec = [dD12dr dD12dz dD12dkr dD12dkz];
        dD22_vec = [dD22dr dD22dz dD22dkr dD22dkz];
        dUdom = dU_vec(1);
        dU_vec   = dU_vec(2:end);

        dU_mat = (1/subsum)*(...
                    dD11_mat*(D22-U) + dD22_mat*(D11-U)...
                    +(dD11_vec - dU_vec).'*(dD22_vec - dU_vec)...
                    +(dD22_vec - dU_vec).'*(dD11_vec - dU_vec)...
                    -2*real(cD12*dD12_mat + dD12_vec'*dD12_vec));

                
        dU_vec = [dUdom dU_vec];
    end



    
case {'cld3x3'}

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
  zroi=zeros(size(Si));                         
  zro=zeros(size(S));
  dLNSids  =dLNnds-2.*omc2./omc2Mom2.*dLNomcds; dSds = sum( Si.*dLNSids ,2);
  dLNSidt  =      -2.*omc2./omc2Mom2.*dLNomcdt; dSdt = sum( Si.*dLNSidt ,2);
  dLNSidom =2*plasma.omega./omc2Mom2;           dSdom= sum( Si.*dLNSidom,2);
  dLNDids  =dLNSids + dLNomcds;                 dDds = sum( Di.*dLNDids ,2);
  dLNDidt  =dLNSidt + dLNomcdt;                 dDdt = sum( Di.*dLNDidt ,2);
  dLNDidom =(3*om2-omc2)./(plasma.omega*omc2Mom2); dDdom= sum( Di.*dLNDidom,2);

  dLNPids  = dLNnds;                             
  dPds     = -sum( Pi.*dLNPids ,2);
  dLNPidt  = zroi;                               
  dPdt     = zro;
  dLNPidom = -2/plasma.omega;                     
  dPdom    = -sum( Pi.*dLNPidom,2);
  
  % second derivatives
  if eval2nd
   omOM = om2 ./omc2Mom2; 
   ocOM = omc2./omc2Mom2;
   dLNSids2 = 2*ocOM.*(2*omOM.*dLNomcds.*dLNomcds -dLNomcds2)+dLNnds2;
   dLNSidst = 2*ocOM.*(2*omOM.*dLNomcds.*dLNomcdt -dLNomcdst);
   dLNSidt2 = 2*ocOM.*(2*omOM.*dLNomcdt.*dLNomcdt -dLNomcdt2);
   dLNDids2 = dLNomcds2 +dLNSids2;
   dLNDidst = dLNomcdst +dLNSidst;
   dLNDidt2 = dLNomcdt2 +dLNSidt2;

   dSids2 = Si.*(dLNSids2 +dLNSids.*dLNSids);   dSds2 = sum( dSids2 ,2);
   dSidst = Si.*(dLNSidst +dLNSids.*dLNSidt);   dSdst = sum( dSidst ,2);
   dSidt2 = Si.*(dLNSidt2 +dLNSidt.*dLNSidt);   dSdt2 = sum( dSidt2 ,2);
   dDids2 = Di.*(dLNDids2 +dLNDids.*dLNDids);   dDds2 = sum( dDids2 ,2);
   dDidst = Di.*(dLNDidst +dLNDids.*dLNDidt);   dDdst = sum( dDidst ,2);
   dDidt2 = Di.*(dLNDidt2 +dLNDidt.*dLNDidt);   dDdt2 = sum( dDidt2 ,2);

   dLNPids2 = dLNnds2;       dLNPidst= zroi;      dLNPidt2 = zroi;
   dPids2   = Pi.*(dLNPids2 +dLNPids.*dLNPids);   dPds2 = -sum( dPids2 ,2);
   dPidst   = Pi.*(dLNPidst +dLNPids.*dLNPidt);   dPdst = -sum( dPidst ,2);
   dPidt2   = Pi.*(dLNPidt2 +dLNPidt.*dLNPidt);   dPdt2 = -sum( dPidt2 ,2);
  end
  
  % --- Tensor elements and derivatives --------------------------------------

    D11 = Nb2+Np2-S;        
    D12 =-Nn.*Nb-D*i;       
    D13 =-Nn.*Np;            
    D22 = Nn2+Np2-S;        
    D23 =-Nb.*Np;             
    D33 = Nn2+Nb2-P;         
    
    cD12 = conj(D12);
    cD13 = conj(D13); 
    cD23 = conj(D23); 
    
    % Dispersion tensor
    DD=[D11 D12 D13 cD12 D22 D23 cD13 cD23 D33];
    % We will assume that the input yv is size [1,4], then reshape works
    DD=reshape(DD,3,3);
    % Convert to double to make sure that double precision LAPACK routines
    % are used.
    [evects,evals]=eig(double(DD));           % Three values
    evs = diag(evals);
    ind = find(abs(evs)==min(abs(evs)));
    U = evs(ind); % The eval nearest zero is the one we want
    pol = evects(:,ind); % get the column with the right evector
    
    % Conversion monitor
    mon2 = abs(det(DD([2 3],[2 3])) + ...
               det(DD([1 3],[1 3])) + ...
               det(DD([1 2],[1 2])) );

    if eval1st
        dD11ds =-dSds;          dD11dt =-dSdt;
        dD12ds =-dDds*i;        dD12dt =-dDdt*i;
        dD13ds = zro;           dD13dt = zro;
        dD22ds =-dSds;          dD22dt =-dSdt;
        dD23ds = zro;           dD23dt = zro;
        dD33ds =-dPds;          dD33dt =-dPdt; 

        dD11dkn= zro;           dD11dkb= Nb*coom*2;     dD11dkp= Np*coom*2;
        dD12dkn=-Nb*coom;       dD12dkb=-Nn*coom;       dD12dkp= zro;
        dD13dkn=-Np*coom;       dD13dkb= zro;           dD13dkp=-Nn*coom;
        dD22dkn= Nn*coom*2;     dD22dkb= zro;           dD22dkp= Np*coom*2;
        dD23dkn= zro;           dD23dkb=-Np*coom;       dD23dkp=-Nb*coom;
        dD33dkn= Nn*coom*2;     dD33dkb= Nb*coom*2;     dD33dkp= zro;

        dD11dkr=dD11dkn.*ener +dD11dkb.*eber +dD11dkp.*eper;
        dD11dkz=dD11dkn.*enez +dD11dkb.*ebez +dD11dkp.*epez;
        dD12dkr=dD12dkn.*ener +dD12dkb.*eber +dD12dkp.*eper;
        dD12dkz=dD12dkn.*enez +dD12dkb.*ebez +dD12dkp.*epez;
        dD13dkr=dD13dkn.*ener +dD13dkb.*eber +dD13dkp.*eper;
        dD13dkz=dD13dkn.*enez +dD13dkb.*ebez +dD13dkp.*epez;
        dD22dkr=dD22dkn.*ener +dD22dkb.*eber +dD22dkp.*eper;
        dD22dkz=dD22dkn.*enez +dD22dkb.*ebez +dD22dkp.*epez;
        dD23dkr=dD23dkn.*ener +dD23dkb.*eber +dD23dkp.*eper;
        dD23dkz=dD23dkn.*enez +dD23dkb.*ebez +dD23dkp.*epez;
        dD33dkr=dD33dkn.*ener +dD33dkb.*eber +dD33dkp.*eper;
        dD33dkz=dD33dkn.*enez +dD33dkb.*ebez +dD33dkp.*epez;

        dD11dom=-2/plasma.omega*(Nb2+Np2)-dSdom;
        dD12dom=-2/plasma.omega*(Nn.*Nb) -dDdom*i;
        dD13dom=-2/plasma.omega*(Nn.*Np);
        dD22dom=-2/plasma.omega*(Nn2+Np2)-dSdom;
        dD23dom=-2/plasma.omega*(Nb.*Np);
        dD33dom=-2/plasma.omega*(Nn2+Nb2)-dPdom;

        DD2=DD-eye(3)*U;

        % Sub-blocks of (dispersion matrix - eigenval*id)
        U1 = det(DD2([2 3],[2 3]));
        U2 = det(DD2([1 3],[1 3]));
        U3 = det(DD2([1 2],[1 2]));
        subsum = U1+U2+U3;

        % Convert to r, z coordinates.
        %  Note that there are extra correction terms due to field
        %  curvature.  They are derivatives of the trasformation elements:
        %  (d/dr)(ebef) etc.
        dNndr = Nr.*denerdr + Nf.*denefdr + Nz.*denezdr;
        dNndz = Nr.*denerdz + Nf.*denefdz + Nz.*denezdz;
        dNbdr = Nr.*deberdr + Nf.*debefdr + Nz.*debezdr;
        dNbdz = Nr.*deberdz + Nf.*debefdz + Nz.*debezdz;
        dNpdr = Nr.*deperdr + Nf.*depefdr + Nz.*depezdr;
        dNpdz = Nr.*deperdz + Nf.*depefdz + Nz.*depezdz;
        
        dD11dr = dD11ds.*dsdr + dD11dt.*dtdr + 2*(Nb.*dNbdr + Np.*dNpdr);
        dD11dz = dD11ds.*dsdz + dD11dt.*dtdz + 2*(Nb.*dNbdz + Np.*dNpdz);
        dD12dr = dD12ds.*dsdr + dD12dt.*dtdr - Nn.*dNbdr - Nb.*dNndr;
        dD12dz = dD12ds.*dsdz + dD12dt.*dtdz - Nn.*dNbdz - Nb.*dNndz;
        dD13dr = dD13ds.*dsdr + dD13dt.*dtdr - dNndr.*Np - Nn.*dNpdr;
        dD13dz = dD13ds.*dsdz + dD13dt.*dtdz - dNndz.*Np - Nn.*dNpdz;
        dD22dr = dD22ds.*dsdr + dD22dt.*dtdr + 2*(Nn.*dNndr + Np.*dNpdr);
        dD22dz = dD22ds.*dsdz + dD22dt.*dtdz + 2*(Nn.*dNndz + Np.*dNpdz);
        dD23dr = dD23ds.*dsdr + dD23dt.*dtdr - dNbdr.*Np - Nb.*dNpdr;
        dD23dz = dD23ds.*dsdz + dD23dt.*dtdz - dNbdz.*Np - Nb.*dNpdz;
        dD33dr = dD33ds.*dsdr + dD33dt.*dtdr + 2*(Nn.*dNndr + Nb.*dNbdr);
        dD33dz = dD33ds.*dsdz + dD33dt.*dtdz + 2*(Nn.*dNndz + Nb.*dNbdz);

        % Gather things together to do all derivatives at the same time
        dD11_vec = [dD11dom dD11dr dD11dz dD11dkr dD11dkz];
        dD12_vec = [dD12dom dD12dr dD12dz dD12dkr dD12dkz];
        dD13_vec = [dD13dom dD13dr dD13dz dD13dkr dD13dkz];
        dD22_vec = [dD22dom dD22dr dD22dz dD22dkr dD22dkz];
        dD23_vec = [dD23dom dD23dr dD23dz dD23dkr dD23dkz];
        dD33_vec = [dD33dom dD33dr dD33dz dD33dkr dD33dkz];

        % Derivatives of the eigenvalue
        %dU_vec   = [dUdom dUdr dUdz dUdkr dUdkz];
        dU_vec = (dD11_vec.*U1 - (D11-U).*2.*real(cD23.*dD23_vec) ...
                +dD22_vec.*U2 - (D22-U).*2.*real(cD13.*dD13_vec) ...
                +dD33_vec.*U3 - (D33-U).*2.*real(cD12.*dD12_vec) )./subsum...
               + 2*real(D12.*D23.*conj(dD13_vec)...
                       +D12.*dD23_vec.*cD13...
                       +dD12_vec.*D23.*cD13)./subsum;
        
    end
    
    % terms needed for second derivative of hamiltonian
    if eval2nd

        dSdr2 = dsdr2.*dSds + dtdr2.*dSdt ...
            + dsdr.*(dsdr.*dSds2 + dtdr.*dSdst) + dtdr.*(dsdr.*dSdst + dtdr.*dSdt2);
        dSdrz = dsdrz.*dSds + dtdrz.*dSdt ...
            + dsdz.*(dsdr.*dSds2 + dtdr.*dSdst) + dtdz.*(dsdr.*dSdst + dtdr.*dSdt2);
        dSdz2 = dsdz2.*dSds + dtdz2.*dSdt ...
            + dsdz.*(dsdz.*dSds2 + dtdz.*dSdst) + dtdz.*(dsdz.*dSdst + dtdz.*dSdt2);
        dDdr2 = dsdr2.*dDds + dtdr2.*dDdt ...
            + dsdr.*(dsdr.*dDds2 + dtdr.*dDdst) + dtdr.*(dsdr.*dDdst + dtdr.*dDdt2);
        dDdrz = dsdrz.*dDds + dtdrz.*dDdt ...
            + dsdz.*(dsdr.*dDds2 + dtdr.*dDdst) + dtdz.*(dsdr.*dDdst + dtdr.*dDdt2);
        dDdz2 = dsdz2.*dDds + dtdz2.*dDdt ...
            + dsdz.*(dsdz.*dDds2 + dtdz.*dDdst) + dtdz.*(dsdz.*dDdst + dtdz.*dDdt2);
        dPdr2 = dsdr2.*dPds + dtdr2.*dPdt ...
            + dsdr.*(dsdr.*dPds2 + dtdr.*dPdst) + dtdr.*(dsdr.*dPdst + dtdr.*dPdt2);
        dPdrz = dsdrz.*dPds + dtdrz.*dPdt ...
            + dsdz.*(dsdr.*dPds2 + dtdr.*dPdst) + dtdz.*(dsdr.*dPdst + dtdr.*dPdt2);
        dPdz2 = dsdz2.*dPds + dtdz2.*dPdt ...
            + dsdz.*(dsdz.*dPds2 + dtdz.*dPdst) + dtdz.*(dsdz.*dPdst + dtdz.*dPdt2);

        dNndr2 = Nr.*denerdr2 + Nf.*denefdr2 + Nz.*denezdr2;
        dNndrz = Nr.*denerdrz + Nf.*denefdrz + Nz.*denezdrz;
        dNndz2 = Nr.*denerdz2 + Nf.*denefdz2 + Nz.*denezdz2;
        dNbdr2 = Nr.*deberdr2 + Nf.*debefdr2 + Nz.*debezdr2;
        dNbdrz = Nr.*deberdrz + Nf.*debefdrz + Nz.*debezdrz;
        dNbdz2 = Nr.*deberdz2 + Nf.*debefdz2 + Nz.*debezdz2;
        dNpdr2 = Nr.*deperdr2 + Nf.*depefdr2 + Nz.*depezdr2;
        dNpdrz = Nr.*deperdrz + Nf.*depefdrz + Nz.*depezdrz;
        dNpdz2 = Nr.*deperdz2 + Nf.*depefdz2 + Nz.*depezdz2;
        
        dD11dr2=2*(dNbdr.*dNbdr + Nb.*dNbdr2 + dNpdr.*dNpdr + Np.*dNpdr2) - dSdr2;
        dD11drz=2*(dNbdr.*dNbdz + Nb.*dNbdrz + dNpdr.*dNpdz + Np.*dNpdrz) - dSdrz;        
        dD11dz2=2*(dNbdz.*dNbdz + Nb.*dNbdz2 + dNpdz.*dNpdz + Np.*dNpdz2) - dSdz2;
        dD12dr2=-(Nn.*dNbdr2 + 2*dNndr.*dNbdr + Nb.*dNndr2)-i*dDdr2;
        dD12drz=-(Nn.*dNbdrz + dNndr.*dNbdz + dNndz.*dNbdr + Nb.*dNndrz)-i*dDdrz;
        dD12dz2=-(Nn.*dNbdz2 + 2*dNndz.*dNbdz + Nb.*dNndz2)-i*dDdz2;
        dD13dr2=-(Nn.*dNpdr2 + 2*dNndr.*dNpdr + Np.*dNndr2);
        dD13drz=-(Nn.*dNpdrz + dNndr.*dNpdz + dNndz.*dNpdr + Np.*dNndrz);
        dD13dz2=-(Nn.*dNpdz2 + 2*dNndz.*dNpdz + Np.*dNndz2);
        dD22dr2=2*(dNndr.*dNndr + Nn.*dNndr2 + dNpdr.*dNpdr + Np.*dNpdr2) - dSdr2;
        dD22drz=2*(dNndr.*dNndz + Nn.*dNndrz + dNpdr.*dNpdz + Np.*dNpdrz) - dSdrz;
        dD22dz2=2*(dNndz.*dNndz + Nn.*dNndz2 + dNpdz.*dNpdz + Np.*dNpdz2) - dSdz2;
        dD23dr2=-(Nb.*dNpdr2 + 2*dNbdr.*dNpdr + Np.*dNbdr2);
        dD23drz=-(Nb.*dNpdrz + dNbdr.*dNpdz + dNbdz.*dNpdr + Np.*dNbdrz);
        dD23dz2=-(Nb.*dNpdz2 + 2*dNbdz.*dNpdz + Np.*dNbdz2);
        dD33dr2=2*(dNndr.*dNndr + Nn.*dNndr2 + dNbdr.*dNbdr + Nb.*dNbdr2) - dPdr2;
        dD33drz=2*(dNndr.*dNndz + Nn.*dNndrz + dNbdr.*dNbdz + Nb.*dNbdrz) - dPdrz;
        dD33dz2=2*(dNndz.*dNndz + Nn.*dNndz2 + dNbdz.*dNbdz + Nb.*dNbdz2) - dPdz2;
        
        dD11drkr = 2*coom*(eber.*dNbdr + Nb.*deberdr + eper.*dNpdr + Np.*deperdr);
        dD11drkz = 2*coom*(ebez.*dNbdr + Nb.*debezdr + epez.*dNpdr + Np.*depezdr);
        dD11dzkr = 2*coom*(eber.*dNbdz + Nb.*deberdz + eper.*dNpdz + Np.*deperdz);
        dD11dzkz = 2*coom*(ebez.*dNbdz + Nb.*debezdz + epez.*dNpdz + Np.*depezdz);
        
        dD12drkr = -coom*(denerdr.*Nb + dNndr.*eber + ener.*dNbdr + Nn.*deberdr);
        dD12drkz = -coom*(denezdr.*Nb + dNndr.*ebez + enez.*dNbdr + Nn.*debezdr);
        dD12dzkr = -coom*(denerdz.*Nb + dNndz.*eber + ener.*dNbdz + Nn.*deberdz);
        dD12dzkz = -coom*(denezdz.*Nb + dNndz.*ebez + enez.*dNbdz + Nn.*debezdz);
        
        dD13drkr = -coom*(denerdr.*Np + dNndr.*eper + ener.*dNpdr + Nn.*deperdr);
        dD13drkz = -coom*(denezdr.*Np + dNndr.*epez + enez.*dNpdr + Nn.*depezdr);
        dD13dzkr = -coom*(denerdz.*Np + dNndz.*eper + ener.*dNpdz + Nn.*deperdz);
        dD13dzkz = -coom*(denezdz.*Np + dNndz.*epez + enez.*dNpdz + Nn.*depezdz);

        dD22drkr = 2*coom*(ener.*dNndr + Nn.*denerdr + eper.*dNpdr + Np.*deperdr);
        dD22drkz = 2*coom*(enez.*dNndr + Nn.*denezdr + epez.*dNpdr + Np.*depezdr);
        dD22dzkr = 2*coom*(ener.*dNndz + Nn.*denerdz + eper.*dNpdz + Np.*deperdz);
        dD22dzkz = 2*coom*(enez.*dNndz + Nn.*denezdz + epez.*dNpdz + Np.*depezdz);

        dD23drkr = -coom*(deperdr.*Nb + dNpdr.*eber + eper.*dNbdr + Np.*deberdr);
        dD23drkz = -coom*(depezdr.*Nb + dNpdr.*ebez + epez.*dNbdr + Np.*debezdr);
        dD23dzkr = -coom*(deperdz.*Nb + dNpdz.*eber + eper.*dNbdz + Np.*deberdz);
        dD23dzkz = -coom*(depezdz.*Nb + dNpdz.*ebez + epez.*dNbdz + Np.*debezdz);

        dD33drkr = 2*coom*(ener.*dNndr + Nn.*denerdr + eber.*dNbdr + Nb.*deberdr);
        dD33drkz = 2*coom*(enez.*dNndr + Nn.*denezdr + ebez.*dNbdr + Nb.*debezdr);
        dD33dzkr = 2*coom*(ener.*dNndz + Nn.*denerdz + eber.*dNbdz + Nb.*deberdz);
        dD33dzkz = 2*coom*(enez.*dNndz + Nn.*denezdz + ebez.*dNbdz + Nb.*debezdz);
        
        dD11dkr2 = 2*coomsq*(eber.^2 + eper.^2);
        dD11dkz2 = 2*coomsq*(ebez.^2 + epez.^2);
        dD11dkrkz= 2*coomsq*(eber.*ebez + eper.*epez);
        dD12dkr2 = -coomsq*(ener.*eber + ener.*eber);
        dD12dkz2 = -coomsq*(ener.*ebez + enez.*eber);
        dD12dkrkz= -coomsq*(enez.*ebez + enez.*ebez);
        dD13dkr2 = -coomsq*(ener.*eper + ener.*eper);
        dD13dkz2 = -coomsq*(ener.*epez + enez.*eper);
        dD13dkrkz= -coomsq*(enez.*epez + enez.*epez);
        dD22dkr2 = 2*coomsq*(ener.^2 + eper.^2);
        dD22dkz2 = 2*coomsq*(enez.^2 + epez.^2);
        dD22dkrkz= 2*coomsq*(ener.*enez + eper.*epez);
        dD23dkr2 = -coomsq*(eper.*eber + eper.*eber);
        dD23dkz2 = -coomsq*(eper.*ebez + epez.*eber);
        dD23dkrkz= -coomsq*(epez.*ebez + epez.*ebez);
        dD33dkr2 = 2*coomsq*(ener.^2 + eber.^2);
        dD33dkz2 = 2*coomsq*(enez.^2 + ebez.^2);
        dD33dkrkz= 2*coomsq*(ener.*enez + eber.*ebez);
        
        % dD11_mat = [ds dt dkr dkz].'*[ds dt dkr dkz];
        dD11_mat = [dD11dr2   dD11drz   dD11drkr  dD11drkz;
                    dD11drz   dD11dz2   dD11dzkr  dD11dzkz;
                    dD11drkr  dD11dzkr  dD11dkr2  dD11dkrkz;
                    dD11drkz  dD11dzkz  dD11dkrkz dD11dkz2 ];
        dD12_mat = [dD12dr2   dD12drz   dD12drkr  dD12drkz;
                    dD12drz   dD12dz2   dD12dzkr  dD12dzkz;
                    dD12drkr  dD12dzkr  dD12dkr2  dD12dkrkz;
                    dD12drkz  dD12dzkz  dD12dkrkz dD12dkz2 ];
        dD13_mat = [dD13dr2   dD13drz   dD13drkr  dD13drkz;
                    dD13drz   dD13dz2   dD13dzkr  dD13dzkz;
                    dD13drkr  dD13dzkr  dD13dkr2  dD13dkrkz;
                    dD13drkz  dD13dzkz  dD13dkrkz dD13dkz2 ];
        dD22_mat = [dD22dr2   dD22drz   dD22drkr  dD22drkz;
                    dD22drz   dD22dz2   dD22dzkr  dD22dzkz;
                    dD22drkr  dD22dzkr  dD22dkr2  dD22dkrkz;
                    dD22drkz  dD22dzkz  dD22dkrkz dD22dkz2 ];
        dD23_mat = [dD23dr2   dD23drz   dD23drkr  dD23drkz;
                    dD23drz   dD23dz2   dD23dzkr  dD23dzkz;
                    dD23drkr  dD23dzkr  dD23dkr2  dD23dkrkz;
                    dD23drkz  dD23dzkz  dD23dkrkz dD23dkz2 ];
        dD33_mat = [dD33dr2   dD33drz   dD33drkr  dD33drkz;
                    dD33drz   dD33dz2   dD33dzkr  dD33dzkz;
                    dD33drkr  dD33dzkr  dD33dkr2  dD33dkrkz;
                    dD33drkz  dD33dzkz  dD33dkrkz dD33dkz2 ];


        % for the second order derivatives, we don't need dom
        dD11_vec = [dD11dr dD11dz dD11dkr dD11dkz];
        dD12_vec = [dD12dr dD12dz dD12dkr dD12dkz];
        dD13_vec = [dD13dr dD13dz dD13dkr dD13dkz];
        dD22_vec = [dD22dr dD22dz dD22dkr dD22dkz];
        dD23_vec = [dD23dr dD23dz dD23dkr dD23dkz];
        dD33_vec = [dD33dr dD33dz dD33dkr dD33dkz];
        dUdom    = dU_vec(1);
        dU_vec   = dU_vec(2:end);

        dU1_vec =  (dD22_vec - dU_vec).*(D33 - U)...
                 + (dD33_vec - dU_vec).*(D22 - U)...
                 -2*real(cD23*dD23_vec);
        dU2_vec =  (dD11_vec - dU_vec).*(D33 - U)...
                 + (dD33_vec - dU_vec).*(D11 - U)...
                 -2*real(cD13*dD13_vec);
        dU3_vec =  (dD22_vec - dU_vec).*(D11 - U)...
                 + (dD11_vec - dU_vec).*(D22 - U)...
                 -2*real(cD12*dD12_vec);

        dB_mat = 2*real( D23.*dD12_vec.'*conj(dD13_vec) ...
                 + D12.*dD23_vec.'*conj(dD13_vec) ...
                 + D12.*D23.*conj(dD13_mat) ...
                 + dD12_vec.'*dD23_vec.*cD13 ...
                 + D12.*dD23_mat.*cD13 ...
                 + D12.*conj(dD13_vec).'*dD23_vec ...
                 + dD12_mat.*D23.*cD13 ...
                 + dD23_vec.'*dD12_vec.*cD13 ...
                 + conj(dD13_vec).'*dD12_vec.*D23 );

        dU_mat = (1/subsum).*( ...
                 dD11_mat*U1 + dD22_mat*U2 + dD33_mat*U3 ...
                 +dU1_vec.'*(dD11_vec-dU_vec) ...
                 +dU2_vec.'*(dD22_vec-dU_vec) ...
                 +dU3_vec.'*(dD33_vec-dU_vec) ...
                 -2*(dD11_vec-dU_vec).'*real(cD23*dD23_vec) ...
                 -2*(dD22_vec-dU_vec).'*real(cD13*dD13_vec) ...
                 -2*(dD33_vec-dU_vec).'*real(cD12*dD12_vec) ...
                 -(D11-U)*2*real(cD23*dD23_mat +dD23_vec'*dD23_vec) ...% ' not .' => conj
                 -(D22-U)*2*real(cD13*dD13_mat +dD13_vec'*dD13_vec) ...
                 -(D33-U)*2*real(cD12*dD12_mat +dD12_vec'*dD12_vec) ...
                 +dB_mat);

        test = dU_mat;
        if nnz( (test - test.') > 1E-10)
            disp(test - test.');
            error('disp_eig:cld3x3:symmetric','Matrix not symmetric');
        end
        
        dU_vec = [dUdom dU_vec];

    end
    

% +++++ Warm plasma model, etc ++++++++++++++++++++++++++++++++++++++++++++++
otherwise
  error(['Unknown model ''' plasma.MODEL(1:6) '''.']);
end;


%
% ===== Variable output dimensions ===========================================
%
switch oper
case {'Mon'}                            % --- Conversion monitors
  mon1=0;
  tmp=real(cat(1,mon1',mon2'));
  out=reshape(tmp,1,numel(tmp));
case {'Dsp','Msw'}                      % --- Dispersion relations
  out = real(U);
case {'Ten'}                            % --- Dispersion tensor
  out = reshape(DD,9,1);
case {'Trj'}                            % --- ODE trajectory
 % dUdom = dU_vec(1);
 %  dUdomM=(1./dUdom)';                 % fn of t, not sigma 
 dUdomM = 1;                            % --- Solve as a fn of sigma, not t
 J = [zeros(2,2) eye(2); -eye(2) zeros(2,2)]; % I think the sign here is different
% dU = [dUdr dUdz dUdkr dUdkz].';
 dU = reshape(dU_vec(2:5),4,1);
 dz = (J*dU);
 if eval2nd
     d2U = dU_mat; 
     ds = reshape(J*d2U*sym_mat,16,1);
     out = real(dUdomM*([dz; ds ]));
 else
     out = real(dUdomM*dz);
 end
   
% case {'Sgn'}
%   out = sign(dUdom); % return the sign of the derivative of time wrt ray parameter
case {'Pol'}
  out = pol;
otherwise
  error(['Unknown oper ''' oper '''.']);
end

function out = cgamma(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Complex gamma function for -10<Re(z)<10, 10<Im(z)<10
%
%  Tested 12 digits accuracy for
%   cgamma([0 .5+.5i -.5+.5i -.5-.5i .5-.5i 1 1+i i -1+i -1 -1-i -i 1-i]')
%
%  A. JAUN, Numerical Analysis, KTH, 100 44 Stockholm, Sweden
%  A.N. KAUFMAN, Lawrence Berkeley Laboratory, Berkeley, CA 94720, USA
%  E.R. TRACY, College of William & Mary, Williamsburg, VA 23187-8795, USA
%
%  Documented under "http://www.nada.kth.se/~jaun"
%
%  (C) Version 7.0,  14-Aug-2006. All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
x=real(z); y=imag(z); n=prod(size(x));
%
% ----- Map to domain of asymptotic expansion
%
for k=1:n
  if (y(k)<0) 
    z(k)=conj(z(k));
  end
end
%
lnCgamma=zeros(size(x));
for k=1:n
  for m=1:floor(9-x)
   lnCgamma(k)=lnCgamma(k)-log(z(k)); z(k)=z(k)+1;
  end
end
%
% ----- First compute log(cgamma) from expansion for 9<Re(z)<10, 0<Im(z)<10
%
zm =1./z;   zm1 =zm; zm2=zm1.*zm1;        % Avoid power function
zm=zm.*zm2; zm3 =zm;
zm=zm.*zm2; zm5 =zm;
zm=zm.*zm2; zm7 =zm;
zm=zm.*zm2; zm9 =zm;
zm=zm.*zm2; zm11=zm;
zm=zm.*zm2; zm13=zm;
zm=zm.*zm2; zm15=zm;
%
lnCgamma=lnCgamma -zm15*(3617./122400);   % Add series in reverse order
lnCgamma=lnCgamma +zm13/156;
lnCgamma=lnCgamma -zm11*(691./360360);
lnCgamma=lnCgamma +zm9 /1188;
lnCgamma=lnCgamma -zm7 /1680;
lnCgamma=lnCgamma +zm5 /1260;
lnCgamma=lnCgamma -zm3 /360;
lnCgamma=lnCgamma +zm1 /12 +(z-0.5).*log(z) -z +0.5*log(2*pi);
%
for k=1:n
%  if (y(k)>0)                             % Don't understand why not <0
  if (y(k)<0)                             % modified by Steve Richardson
    lnCgamma(k)=conj(lnCgamma(k));
  end
end
%
% ----- Take exponential
%
out = exp(lnCgamma);

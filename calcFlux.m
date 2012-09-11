function plasma = calcFlux(plasma)

    %% Calculate the flux surface coordinates ('mesh' and 'boundary')
    s=0.0001:(.999/(plasma.NS-1)):1.;                 %  radial mesh
    theta=0:2*pi/plasma.NT:2*pi;                      %  poloidal mesh
    plasma.r = zeros(plasma.NS,plasma.NT+1);
    plasma.z = zeros(plasma.NS,plasma.NT+1);
    for kk=1:plasma.NS
        si=s(kk).*ones(size(theta));             %  flux surfaces
        [rho,plasma.r(kk,:),plasma.z(kk,:)]=mapFlux(si,theta); %    coordinates
    end;

    plasma.Rho = sqrt((plasma.r-plasma.r0).^2 +plasma.z.^2);
    plasma.Theta = atan2(-plasma.z,(plasma.r-plasma.r0));

    

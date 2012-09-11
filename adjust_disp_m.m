function y_ant = adjust_disp_m(y,m)
    global plasma %cnst
    
    % ----- Adjusted position & wave vector in cylindrical coordinates -----------
    r=y(1); z=y(2); 
    kr=y(3); kz=y(4); % starting wave vector
    rho=sqrt((r-plasma.r0).^2+z.^2);
    theta=atan2(plasma.r0-r,-z)+pi/2;

    krho_old=-sqrt(kr^2 + kz^2 -(m/rho)^2); % assume wave starts out going in
    %krho=+sqrt(kr^2 + kz^2 -(m/rho)^2); 

    % Want to adjust kr and kz so that, for given m, k_phi, the initial
    % condtions are on the dispersion surface.
    % m = function argument = rho*k_theta
    
    %kf = plasma.kant(2);

    %y0=[r,z,kr,(kr*sin(theta)+(m/rho))/cos(theta)];
    
    % kz = (kr*sin(theta)+(m/rho))/cos(theta)
%     f = @(krho) dispersion([r,z,...
%         krho*cos(theta)-(m/rho)*sin(theta),...
%         -(krho*sin(theta)+(m/rho)*cos(theta))],'Dsp');
%     krho = fzero(f,krho);
    
    f = @(krho) disp_eig([r,z,...
        krho*cos(theta)-(m/rho)*sin(theta),...
        -(krho*sin(theta)+(m/rho)*cos(theta)), reshape(eye(4),1,16)].','Dsp');
%     f = @(krho) disp_eig([r,z,...
%         krho*cos(theta)-(m/rho)*sin(theta),...
%         -(krho*sin(theta)+(m/rho)*cos(theta)), zeros(1,16)].','Dsp');
    krho_new = fzero(f,krho_old);
    
    %kz = (kr*sin(theta)+(m/rho))/cos(theta);
    
    y_ant = [r,z,krho_new*cos(theta)-(m/rho)*sin(theta),...
             -(krho_new*sin(theta)+(m/rho)*cos(theta))].';
    
    % ----- Wave vector and refraction index -------------------------------------
%     kn=kr.*ener +kf.*enef +kz.*enez;
%     kb=kr.*eber +kf.*ebef +kz.*ebez;
%     kp=kr.*eper +kf.*epef +kz.*epez;
%     T =[ener enef enez;...
%         eber ebef ebez;...
%         eper epef epez];
%     k_cart = T^(-1) * [kn;kb;kp];
end


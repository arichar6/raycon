function plasma = plotFlux(plasma) 
    grey = [0.6200    0.6200    0.6200];
    
    if isempty(plasma.r) || isempty(plasma.z)
        plasma = calcFlux(plasma);
    end
    
    % Plot the flux surface coordinates ('mesh' and 'boundary')
    h=plot(plasma.r(:,:),plasma.z(:,:)); % Plot radial lines
    set(h,'Color',grey); 
    hold on
    % Plot flux surfaces
    h=plot(plasma.r(:,:).',plasma.z(:,:).');
    set(h,'Color',grey); 
    h=plot(plasma.r(end,:),plasma.z(end,:)); 
    % Make the last fulx surface a thicker line
    set(h,'LineWidth',2,'Color',grey); 
    hold off
    xlabel('major radius  R [m]');  zlabel('vertical direction  Z [m]');
    %axis('equal');axis('tight');
    
end
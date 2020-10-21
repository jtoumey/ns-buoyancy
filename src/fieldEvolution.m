%--------------------------------------------------------------------------
% ME5311, Spring 2020
% Term Project, Part 3
% Julian M. Toumey
% 08 May 2020
%--------------------------------------------------------------------------
clc

%%
dirTheta = dir('**/theta_t0*.dat');
dirUvel = dir('**/uVel_t*.dat');
dirVvel = dir('**/vVel_t*.dat');

xNodeLoc = load("nodeLocation_x.dat");
yNodeLoc = load("nodeLocation_y.dat");

% Interpolate node coordinates to cell centers
xNodeLoc = xNodeLoc(1:end-1, 1:end-1);
dx = xNodeLoc(2, 1) - xNodeLoc(1, 1);
xNodeLoc = xNodeLoc + dx/2;

yNodeLoc = yNodeLoc(1:end-1, 1:end-1);
dy = yNodeLoc(1, 2) - yNodeLoc(1, 1);
yNodeLoc = yNodeLoc + dy/2;

% Animation details
h = figure;
filename = "buoyancyAnim.gif";
colormap cool; % hot/winter/summer

for nn = 1:numel(dirTheta)
    % Theta
    th = load(dirTheta(nn).name);
    subplot(2, 1, 1);
     
    h1 = pcolor(xNodeLoc, yNodeLoc, th);
    set(h1, 'EdgeColor', 'none');
    colorbar;
    caxis([0 1]);

    % Velocity
    u0 = load(dirUvel(nn).name);
    u0 = u0(2:end, :) - u0(1:end-1, :);
    v0 = load(dirVvel(nn).name);
    v0 = v0(:, 2:end) - v0(:, 1:end-1);

    subplot(2,1,2);
    quiver(xNodeLoc, yNodeLoc, u0, v0);
    
    pause(0.1);
    
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if nn == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
    
end
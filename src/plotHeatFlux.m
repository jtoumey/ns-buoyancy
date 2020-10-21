%--------------------------------------------------------------------------
% ME5311, Spring 2020
% Term Project, Part 3
% Julian M. Toumey
% 08 May 2020
%--------------------------------------------------------------------------
clc

%%
dirTheta = dir('**/theta_t0*.dat');

xNodeLoc = load("nodeLocation_x.dat");
yNodeLoc = load("nodeLocation_y.dat");

% Interpolate node coordinates to cell centers
xNodeLoc = xNodeLoc(1:end-1, 1:end-1);
dx = xNodeLoc(2, 1) - xNodeLoc(1, 1);
xNodeLoc = xNodeLoc + dx/2;

yNodeLoc = yNodeLoc(1:end-1, 1:end-1);
dy = yNodeLoc(1, 2) - yNodeLoc(1, 1);
yNodeLoc = yNodeLoc + dy/2;

dx = xNodeLoc(2, 1) - xNodeLoc(1, 1);
dy = yNodeLoc(1, 2) - yNodeLoc(1, 1);

thWbc = 1.0;
thEbc = 0.0;

gradTw = zeros(1, numel(dirTheta));
gradTe = zeros(1, numel(dirTheta));

for nn = 1:numel(dirTheta)
    % Load each Theta field
    th = load(dirTheta(nn).name);
    
    % Compute each wall heat flux
    thetaWest = th(1, :);
    thetaEast = th(end, :);
    gradTw(nn) = -mean((thetaWest - thWbc)/dx);
    gradTe(nn) = -mean((thEbc - thetaEast)/dy);
    
end

% Visualize results
plot(gradTw, '--', 'LineWidth', 2)
hold on;
plot(gradTe, '-.', 'LineWidth', 2)

legend('West Wall', 'East Wall', 'Interpreter', 'latex', 'fontsize', 12);
title('Average Non-Dim Wall Heat Flux', 'Interpreter', 'latex', 'fontsize', 14);
xlabel('Time-step', 'Interpreter', 'latex', 'fontsize', 12);
ylabel('Non-Dimensional Wall Heat Flux', 'Interpreter', 'latex', 'fontsize', 12);


%% Save tha fig
fig = gcf;
% fig.PaperPositionMode = 'auto';
orient(fig, 'landscape');
print('complexHeatFlux.pdf','-dpdf','-fillpage')
% saveas(gcf, 'heatFlux.pdf')


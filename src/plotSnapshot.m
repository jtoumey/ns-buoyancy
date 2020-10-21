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

colormap cool;

% Pick some times to plot
dTs = floor(numel(dirTheta)/3);
t0 = 1;
t1 = dTs;
t2 = dTs*2;
t3 = numel(dirTheta);
timeSteps = [t0 t1 t2 t3];

%% Time 1: t = 0
nn = timeSteps(1);

th0 = load(dirTheta(nn).name);
subplot(3, 4, 1);
h1 = pcolor(xNodeLoc, yNodeLoc, th0);
set(h1, 'EdgeColor', 'none');
ylabel("$\Theta$", 'Interpreter', 'latex');

% Parse the file to get the current time
tt = dirTheta(nn).name;
tt = strsplit(tt, '_t');
tt = char(tt(2));
tt = strsplit(tt, '.d');
tt = char(tt(1));
titleString = strcat('Time $t = ', tt, '$ s');
title(titleString, 'Interpreter', 'latex');

% u-Velocity
u0 = load(dirUvel(nn).name);
u0 = u0(2:end, :) - u0(1:end-1, :);
subplot(3, 4, 5)
h2 = pcolor(xNodeLoc, yNodeLoc, u0);
set(h2, 'EdgeColor', 'none');
ylabel({"Non-dimensional", "$u$ Velocity"}, 'Interpreter', 'latex');

% v-Velocity
v0 = load(dirVvel(nn).name);
v0 = v0(:, 2:end) - v0(:, 1:end-1);
subplot(3, 4, 9)
h3 = pcolor(xNodeLoc, yNodeLoc, v0);
set(h3, 'EdgeColor', 'none');
ylabel({"Non-dimensional", "$v$ Velocity"}, 'Interpreter', 'latex');

%% Time 2: t = 
nn = timeSteps(2);

th0 = load(dirTheta(nn).name);
subplot(3, 4, 2);
h1 = pcolor(xNodeLoc, yNodeLoc, th0);
set(h1, 'EdgeColor', 'none');

% Parse the file to get the current time
tt = dirTheta(nn).name;
tt = strsplit(tt, '_t');
tt = char(tt(2));
tt = strsplit(tt, '.d');
tt = char(tt(1));
titleString = strcat('Time $t = ', tt, '$ s');
title(titleString, 'Interpreter', 'latex');


% u-Velocity
u0 = load(dirUvel(nn).name);
u0 = u0(2:end, :) - u0(1:end-1, :);
subplot(3, 4, 6)
h2 = pcolor(xNodeLoc, yNodeLoc, u0);
set(h2, 'EdgeColor', 'none');

% v-Velocity
v0 = load(dirVvel(nn).name);
v0 = v0(:, 2:end) - v0(:, 1:end-1);
subplot(3, 4, 10)
h3 = pcolor(xNodeLoc, yNodeLoc, v0);
set(h3, 'EdgeColor', 'none');


%% Time 3: t = 
nn = timeSteps(3);

th0 = load(dirTheta(nn).name);
subplot(3, 4, 3);
h1 = pcolor(xNodeLoc, yNodeLoc, th0);
set(h1, 'EdgeColor', 'none');

% Parse the file to get the current time
tt = dirTheta(nn).name;
tt = strsplit(tt, '_t');
tt = char(tt(2));
tt = strsplit(tt, '.d');
tt = char(tt(1));
titleString = strcat('Time $t = ', tt, '$ s');
title(titleString, 'Interpreter', 'latex');

% u-Velocity
u0 = load(dirUvel(nn).name);
u0 = u0(2:end, :) - u0(1:end-1, :);
subplot(3, 4, 7)
h2 = pcolor(xNodeLoc, yNodeLoc, u0);
set(h2, 'EdgeColor', 'none');

% v-Velocity
v0 = load(dirVvel(nn).name);
v0 = v0(:, 2:end) - v0(:, 1:end-1);
subplot(3, 4, 11)
h3 = pcolor(xNodeLoc, yNodeLoc, v0);
set(h3, 'EdgeColor', 'none');

%% Time 4: t = 
nn = timeSteps(4);

th0 = load(dirTheta(nn).name);
subplot(3, 4, 4);
h1 = pcolor(xNodeLoc, yNodeLoc, th0);
set(h1, 'EdgeColor', 'none');

% Parse the file to get the current time
tt = dirTheta(nn).name;
tt = strsplit(tt, '_t');
tt = char(tt(2));
tt = strsplit(tt, '.d');
tt = char(tt(1));
titleString = strcat('Time $t = ', tt, '$ s');
title(titleString, 'Interpreter', 'latex');

hp4 = get(subplot(3, 4, 4),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(4)])
caxis([0 1]);

% u-Velocity
u0 = load(dirUvel(nn).name);
u0 = u0(2:end, :) - u0(1:end-1, :);
subplot(3, 4, 8)
h2 = pcolor(xNodeLoc, yNodeLoc, u0);
set(h2, 'EdgeColor', 'none');

hp4 = get(subplot(3, 4, 8),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(4)])

% v-Velocity
v0 = load(dirVvel(nn).name);
v0 = v0(:, 2:end) - v0(:, 1:end-1);
subplot(3, 4, 12)
h3 = pcolor(xNodeLoc, yNodeLoc, v0);
set(h3, 'EdgeColor', 'none');

hp4 = get(subplot(3, 4, 12),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(4)])

%% Save the fig
fig = gcf;
fig.PaperPositionMode = 'auto';
print('complexIc.pdf','-dpdf','-fillpage')
% saveas(gcf, 'complexIc.pdf')

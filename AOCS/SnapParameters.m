%% B-Dot High Pass algorithm
% Description: Simulates CubeSat dynamics in low Earth orbit using a modified B-dot control algorithm with a high-pass filter.
% Author:
% Created:

%% Orbit Parameters
REarth = 6371e3; %[m]
h = 425e3; %[m]
mu = 3.986e14; %[m^3/s^2]

alt = REarth + h; %[m]
T = 2 * pi * sqrt((alt ^ 3) / mu); %[s]
w = 2 * pi / T; %[rad/s]

ecc = 0;
inc = 97.787; %[deg]
aop = 144.8873; %[deg]
raan = 59.4276; %[deg]
numOrbits = 3; %[-] Simulate up to 3 orbits

%% CubeSat Parameters
length = 0.1; %[m]
width = 0.1; %[m]
height = 0.2; %[m]
weight = 2.6; %[kg]

% Moments of Inertia (kg*m^2)
I = diag([0.0046, 0.0046, 0.00145]);
invI = inv(I);

yaw0 = 0; %[deg]
pitch0 = 0; %[deg]
roll0 = 0; %[deg]

w0 = [10, 10, 10];
Q0 = [1, 0, 0, 0];

% Magnetorquer Parameters
bDotGain = [1e-4; 1e-4; 1e-4]; % Significantly reduced gain values for stability

% Simulation Parameters
startDateTime = datetime(2019, 5, 19, 12, 0, 0); % Start date and time
dt = 0.1; % Time step based on sample rate
totalTime = numOrbits * T; % Total simulation time in seconds for specified orbits
numSteps = ceil(totalTime / dt) + 1; % Total number of time steps, ensure integer
timeNumeric = (0:numSteps-1) * dt; % Numeric time array for indexing
timeDatetime = startDateTime + seconds(timeNumeric); % Datetime array for logging and plotting

% Detumble condition
detumbleThreshold = 2.5 * pi / 180; % Angular velocity threshold in rad/s
detumbleAchieved = false; % Flag to indicate if detumbling is achieved

% Magnetometer
Imax = 0.2143;

%% Computing the magnetic field
% Latitude and Longitude
lat0 = 40.33175117884342;
long0 = -3.7667574322954365;

K_Matrix = [0.5,0,0;0,0.5,0;0,0,0.2];


%% Magnetometer parameters

linearity_mmeter = 0.005;
bias = [28.8675 28.8675 28.8675]; %[nT]
noise = 0; %[nT]
noise_variance = 30; %[nT]
missalingment_x_mmeter = 0.5*pi/180; %[rad]
missalingment_y_mmeter = 0.5*pi/180; %[rad]
missalingment_z_mmeter = 0.5*pi/180; %[rad]
rotm_mmeter = eul2rotm([missalingment_x_mmeter missalingment_y_mmeter missalingment_z_mmeter]); %Missalingment rotation matrix
% Non-linearities lookup table calculation using sinusoidal approach
lookup_breakpoints_mmeter = linspace(-200,200,50)*10^3; %[nT]
max_breakpoint_mmeter = max(lookup_breakpoints_mmeter);
max_deviation_mmeter = max_breakpoint_mmeter*linearity_mmeter;
lookup_tabledata_mmeter = lookup_breakpoints_mmeter+max_deviation_mmeter*sin(pi*lookup_breakpoints_mmeter/max_breakpoint_mmeter); %[nT]

%% Magnetorquer

% Maximum magnetic torque that can be produced by the torquers and the coil
boira_max_magnetic_torque = [0.2 0.2 0.12]; % (A m^2) (body frame)
boira_magnetic_gains = [2.3 2.3 2.1]; % Magnetic gains of the manetorquers (A*m²/A)

%Maximum currents to circulate through the x, y and z magnetorquers (A)
boira_Imax = boira_max_magnetic_torque./boira_magnetic_gains; 

% Linearity of the magnetorquers in the relation intensity-magnetic moment
boira_linearity = [2.5 2.5 2.5]/100; 

missalingment_x_mtorquer = 1*pi/180; %[rad]
missalingment_y_mtorquer = 1*pi/180; %[rad]
missalingment_z_mtorquer = 1*pi/180; %[rad]
rotm_mtorquer = eul2rotm([missalingment_x_mtorquer missalingment_y_mtorquer missalingment_z_mtorquer]); %Missalingment rotation matrix

% Non-linearities lookup table calculation using sinusoidal approach
linearity_mtorquer = boira_linearity(1);
lookup_breakpoints_mtorquer = linspace(-boira_Imax(1)*2,boira_Imax(1)*2,55);
max_breakpoint_mtorquer = max(lookup_breakpoints_mtorquer);
max_deviation_mtorquer = max_breakpoint_mtorquer*linearity_mtorquer;
lookup_tabledata_mtorquer = lookup_breakpoints_mtorquer+max_deviation_mtorquer*sin(pi*lookup_breakpoints_mtorquer/max_breakpoint_mtorquer); 

linearity_mcoil = boira_linearity(3);
lookup_breakpoints_mcoil = linspace(-boira_Imax(3)*2,boira_Imax(3)*2,60);
max_breakpoint_mcoil = max(lookup_breakpoints_mcoil);
max_deviation_mcoil = max_breakpoint_mcoil*linearity_mcoil;  
lookup_tabledata_mcoil = lookup_breakpoints_mcoil+max_deviation_mcoil*sin(pi*lookup_breakpoints_mcoil/max_breakpoint_mcoil);


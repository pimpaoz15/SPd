%% B-Dot High Pass algorithm
% Description: Simulates CubeSat dynamics in low Earth orbit using a modified B-dot control algorithm with a high-pass filter.
% Author:
% Created:

clear;
clc;

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

% High-pass filter cutoff frequency
f_c = 1; % Example value, adjust as needed

% Detumble condition
detumbleThreshold = 2.5 * pi / 180; % Angular velocity threshold in rad/s
detumbleAchieved = false; % Flag to indicate if detumbling is achieved

% Magnetometer
Imax = 0.2143;

%% Computing the magnetic field
% Latitude and Longitude
lat0 = 40.33175117884342;
long0 = -3.7667574322954365;
%% B-Dot High Pass Algorithm
% Description: Simulates CubeSat dynamics in low Earth orbit using a modified B-dot control algorithm with a high-pass filter.
% Author:
% Created:

clear;
clc;

%% Orbit Parameters
REarth = 6371e3; %[m]
h = 425e3; %[m]
mu = 3.986e14; %[m^3/s^2]

altitude = REarth + h; %[m]
T = 2 * pi * sqrt((altitude ^ 3) / mu); %[s]
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

% Moment of Inertia (kg*m^2)
I = diag([0.0046, 0.0046, 0.00145]);

% Magnetorquer Parameters
bDotGain = [1e-4; 1e-4; 1e-4]; % Significantly reduced gain values for stability

% Simulation Parameters
dt = 1; % Reduced time step in seconds
totalTime = numOrbits * T; % Total simulation time in seconds for specified orbits
numSteps = ceil(totalTime / dt) + 1; % Total number of time steps, ensure integer
timeNumeric = (0:numSteps-1) * dt; % Numeric time array for indexing

% High-pass filter cutoff frequency
f_c = 1; % Example value, adjust as needed

% Detumble condition
detumbleThreshold = 2.5 * pi / 180; % Angular velocity threshold in rad/s
detumbleAchieved = false; % Flag to indicate if detumbling is achieved

%% Computing the Magnetic Field
% Latitude and Longitude
latitude = 40.33175117884342;
longitude = -3.7667574322954365;

B = igrfmagm(altitude, latitude, longitude, 2024);
B = B(:); % Ensure B is a column vector

% Define constants
deg2rad = pi / 180;
rad2deg = 180 / pi;
max_orbits = 4; % Maximum orbits to complete detumbling
angular_velocity_threshold = 2.5 * deg2rad; % Threshold in rad/s

% Initial conditions
initial_angular_velocity = [3, 0, 0] * deg2rad; % Initial angular velocity in rad/s

% Initialize orbit and state variables
current_angular_velocity = initial_angular_velocity(:); % Ensure it's a column vector
prev_B = zeros(3, 1); % Initialize previous magnetic field

% Arrays to store angular rates for plotting
angular_velocity_history = zeros(numSteps, 3);
time_history = zeros(numSteps, 1);

% Detumbling loop
for step = 1:numSteps
    t = (step - 1) * dt;
    % Propagate orbit
    current_time = t; % Assuming time is in seconds since epoch
    [positions, velocities] = propagateOrbit(current_time, altitude, ecc, inc, raan, aop, 0);
    position = positions(end, :); % Get the latest position
    velocity = velocities(end, :); % Get the latest velocity

    % Convert ECEF to Geodetic coordinates
    TOL = 1e-12;
    [lat, lon, alt] = ecef2geod(positions(1), positions(2), positions(3), TOL);

    % Calculate magnetic field vector and its derivative (B-dot)
    B = igrfmagm(alt, lat, lon, 2024); % Ensure B is a column vector
    B = B(:); % Ensure B is a column vector
    B_dot = (B - prev_B) / dt;
    prev_B = B;

    % B-dot control law
    K = -100; % Reduced control gain for stability
    magnetic_moment = K * B_dot(:); % Ensure magnetic_moment is a column vector

    % Update angular velocity (simple model, assuming magnetic torquers)
    torque = cross(magnetic_moment, B);
    angular_acceleration = inv(I) * torque;
    current_angular_velocity = current_angular_velocity + angular_acceleration * dt;

    % Store angular velocity for plotting
    angular_velocity_history(step, :) = current_angular_velocity(:)';
    time_history(step) = t;

    % Display step counter
    disp(['Step: ', num2str(t/dt), ', Time: ', num2str(t), ' s']);

    % Check if detumbling criteria are met
    if norm(current_angular_velocity) < angular_velocity_threshold
        disp('Detumbling successful');
        break;
    end
end

if norm(current_angular_velocity) >= angular_velocity_threshold
    disp('Detumbling not completed within 4 orbits');
end

% Plot angular rates
figure;
subplot(3, 1, 1);
plot(time_history, angular_velocity_history(:, 1) * rad2deg);
xlabel('Time (s)');
ylabel('Angular Rate X (deg/s)');
title('Angular Rate X');

subplot(3, 1, 2);
plot(time_history, angular_velocity_history(:, 2) * rad2deg);
xlabel('Time (s)');
ylabel('Angular Rate Y (deg/s)');
title('Angular Rate Y');

subplot(3, 1, 3);
plot(time_history, angular_velocity_history(:, 3) * rad2deg);
xlabel('Time (s)');
ylabel('Angular Rate Z (deg/s)');
title('Angular Rate Z');


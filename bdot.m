%% B-Dot High Pass Algorithm
% Description: Simulates CubeSat dynamics in low Earth orbit using a modified B-dot control algorithm with a high-pass filter.
% Author:
% Created:

clear;
clc;
close all;

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
numOrbits = 4; %[-] Simulate up to 4 orbits

%% CubeSat Parameters
length = 0.1; %[m]
width = 0.1; %[m]
height = 0.2; %[m]
weight = 2.6; %[kg]

% Moment of Inertia (kg*m^2)
I = diag([0.0046, 0.0046, 0.00145]);
invI = inv(I);

% Magnetorquer Parameters
bDotGain = [1e-4; 1e-4; 1e-4]; % Adjusted gain values for stability

% Simulation Parameters
dt = 0.1; % Time step in seconds
totalTime = numOrbits * T; % Total simulation time in seconds for specified orbits
numSteps = ceil(totalTime / dt) + 1; % Total number of time steps, ensure integer
timeNumeric = (0:numSteps-1) * dt; % Numeric time array for indexing

% High-pass filter cutoff frequency
f_c = 1; % Example value, adjust as needed

% Detumble condition
detumbleThreshold = 2.5 * pi / 180; % Angular velocity threshold in rad/s
detumbleAchieved = false; % Flag to indicate if detumbling is achieved

%% Initialize Variables
latitude = 40.33175117884342;
longitude = -3.7667574322954365;
timeDatetime = datetime(2024, 1, 1, 0, 0, 0) + seconds(timeNumeric);

B = zeros(3, numSteps);
Bdot = zeros(3, numSteps);
controlTorques = zeros(3, numSteps);
angularRates = zeros(3, numSteps);
attitudes = zeros(3, numSteps);

% Initial conditions
initial_angular_velocity = [3, 0, 0] * (pi / 180); % Initial angular velocity in rad/s
angularRates(:, 1) = initial_angular_velocity;

% Initial magnetic field
Bprev = igrfmagm(altitude, latitude, longitude, 2024);
Bprev = Bprev(:); % Ensure Bprev is a column vector
B(:, 1) = Bprev;

% Initial position and velocity
[r0, v0] = propagateOrbit(timeDatetime(1), altitude, ecc, inc, raan, aop, 0);

%% Detumbling Loop
for k = 2:numSteps
    % Propagate orbit to get current position and velocity
    [positions, velocities] = propagateOrbit(timeDatetime(1:k), altitude, ecc, inc, raan, aop, 0);
    
    % Get the current position and velocity
    position = positions(:, end);
    velocity = velocities(:, end);
    
    % Convert ECEF to Geodetic coordinates
    TOL = 1e-12;
    [lat, lon, alt] = ecef2geod(position(1), position(2), position(3), TOL);
    
    % Update magnetic field using IGRF model
    currentDatetime = timeDatetime(k);
    B_current = igrfmagm(alt, lat, lon, 2024);
    B_current = B_current(:); % Ensure B_current is a column vector
    B(:, k) = B_current;
    
    % High-pass filter for B-dot calculation
    Bdot(:, k) = (B_current - Bprev) / dt; % Derivative of B
    Bprev = B_current;

    % B-dot Control Algorithm active only during the first 3 orbits or until detumbling is achieved
    currentOrbit = floor(timeNumeric(k) / T) + 1;
    if currentOrbit <= 3 && ~detumbleAchieved
        controlTorques(:, k) = -bDotGain .* Bdot(:, k);

        % Check if detumbling is achieved
        if norm(angularRates(:, k-1)) < detumbleThreshold
            fprintf('Detumbling achieved at step %d (Time: %s, Orbit: %d)\n', k, datestr(currentDatetime), currentOrbit);
            detumbleAchieved = true;
        end
    else
        controlTorques(:, k) = [0; 0; 0]; % No control torques after detumbling phase
    end

    % Compute angular rate derivatives using Euler's equations
    dOmega = eulerEquations(I, angularRates(:, k-1), controlTorques(:, k));

    % Update angular rates using Runge-Kutta integration
    k1 = dOmega;
    k2 = eulerEquations(I, angularRates(:, k-1) + 0.5 * k1 * dt, controlTorques(:, k));
    k3 = eulerEquations(I, angularRates(:, k-1) + 0.5 * k2 * dt, controlTorques(:, k));
    k4 = eulerEquations(I, angularRates(:, k-1) + k3 * dt, controlTorques(:, k));
    angularRates(:, k) = angularRates(:, k-1) + (1/6) * (k1 + 2*k2 + 2*k3 + k4) * dt;

    % Check for numerical stability
    if any(isnan(angularRates(:, k))) || any(isinf(angularRates(:, k)))
        warning('Numerical instability detected at time step %d. Adjusting time step or parameters might be necessary.', k);
        break;
    end

    % Update attitude (Euler angles)
    attitudes(:, k) = attitudes(:, k-1) + angularRates(:, k) * dt;

    % Print progress to console every 100 steps
    if mod(k, 100) == 0
        fprintf('Step %d/%d, Time: %s, Orbit: %d\n', k, numSteps, datestr(currentDatetime), currentOrbit);
    end
end

% Check if detumbling not completed within 4 orbits
if ~detumbleAchieved
    fprintf('Detumbling not completed within 4 orbits.\n');
end

% Plot angular rates
figure;
subplot(3, 1, 1);
plot(timeNumeric, angularRates(1, :) * (180 / pi));
xlabel('Time (s)');
ylabel('Angular Rate X (deg/s)');
title('Angular Rate X');

subplot(3, 1, 2);
plot(timeNumeric, angularRates(2, :) * (180 / pi));
xlabel('Time (s)');
ylabel('Angular Rate Y (deg/s)');
title('Angular Rate Y');

subplot(3, 1, 3);
plot(timeNumeric, angularRates(3, :) * (180 / pi));
xlabel('Time (s)');
ylabel('Angular Rate Z (deg/s)');
title('Angular Rate Z');

%% Euler's Equations Function
function dOmega = eulerEquations(I, omega, torque)
    dOmega = inv(I) * (torque - cross(omega, I * omega));
end

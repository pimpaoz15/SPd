%% Modified B-Dot High Pass Algorithm with IGRF Magnetic Field
% Description: Simulates CubeSat dynamics in low Earth orbit using a modified B-dot control algorithm with a high-pass filter and IGRF magnetic field.
% Author:
% Created:

clear; clc; close all;

%% Orbit Parameters
REarth = 6371e3; %[m]
h = 425e3; %[m]
mu = 3.986e14; %[m^3/s^2]

alt = REarth + h; %[m]
T = 2 * pi * sqrt((alt ^ 3) / mu); %[s]
w_orbit = 2 * pi / T; %[rad/s]

a = alt; % Semi-major axis
ecc = 0; % Eccentricity
inc = 97.787; %[deg] Inclination
raan = 59.4276; %[deg] Right ascension of ascending node
aop = 144.8873; %[deg] Argument of periapsis
nu = 0; %[deg] True anomaly
numOrbits = 4; %[-] Simulate up to 4 orbits

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

w0 = [10, 10, 10]; % Initial angular velocity [deg/s]
w0 = w0 * pi / 180; % Convert to rad/s
Q0 = [1, 0, 0, 0]; % Initial quaternion

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

% Magnetic Field Parameters
Imax = 0.2143; % Maximum current for magnetorquers

%% Preallocate Arrays
wHistory = zeros(3, numSteps);
QHistory = zeros(4, numSteps);
controlTorqueHistory = zeros(3, numSteps);
bFieldHistory = zeros(3, numSteps);
rHistory = zeros(3, numSteps);
vHistory = zeros(3, numSteps);

%% Initial Conditions
wHistory(:, 1) = w0;
QHistory(:, 1) = Q0;

%% Calculate Initial Position and Velocity
[r, v] = keplerian_to_eci(a, ecc, inc, raan, aop, nu, mu);
rHistory(:, 1) = r;
vHistory(:, 1) = v;

%% Simulation Loop
for k = 1:numSteps-1
    % Propagate position and velocity
    [r, v] = propagateOrbit([0 dt], rHistory(:, k), vHistory(:, k));
    rHistory(:, k+1) = r(:, end);
    vHistory(:, k+1) = v(:, end);
    
    % Convert position to latitude, longitude, and altitude
    [lat, lon, alt] = ecef2lla(rHistory(:, k+1)' * 1e-3); % Convert to km for ecef2lla
    
    % Compute magnetic field using IGRF model
    dyear = decyear(startDateTime.Year, startDateTime.Month, startDateTime.Day + timeNumeric(k)/86400);
    [XYZ, H, D, I, F] = igrfmagm(alt*1e3, lat, lon, dyear);
    bField = XYZ' * 1e-9; % Convert from nT to T (Tesla)
    bFieldHistory(:, k) = bField;
    
    % Compute control torque using modified B-dot
    if k > 1
        bDot = (bField - bFieldHistory(:, k-1)) / dt;
    else
        bDot = [0; 0; 0];
    end
    bHat = bField / norm(bField);
    controlTorque = -bDotGain .* cross(bHat, bDot);
    
    % Store control torque
    controlTorqueHistory(:, k) = controlTorque;

    % Update angular velocity
    w = wHistory(:, k);
    wDot = invI * (controlTorque - cross(w, I*w));
    wHistory(:, k+1) = w + wDot * dt;

    % Update quaternion
    Q = QHistory(:, k);
    QDot = 0.5 * [-Q(2), -Q(3), -Q(4); Q(1), -Q(4), Q(3); Q(4), Q(1), -Q(2); -Q(3), Q(2), Q(1)] * w;
    QHistory(:, k+1) = Q + QDot * dt;
    QHistory(:, k+1) = QHistory(:, k+1) / norm(QHistory(:, k+1)); % Normalize quaternion

    fprintf("here");

    % Check detumble condition
    if norm(wHistory(:, k+1)) < detumbleThreshold
        detumbleAchieved = true;
        fprintf('Detumbling achieved at time step %d (%.2f seconds)\n', k, timeNumeric(k));
        break;
    end
end

%% Results
% Plot Angular Velocity
figure;
plot(timeNumeric(1:k), wHistory(:, 1:k) * 180 / pi); % Convert to deg/s for plotting
xlabel('Time [s]');
ylabel('Angular Velocity [deg/s]');
legend('wx', 'wy', 'wz');
title('CubeSat Angular Velocity with Modified B-dot Control');

% Plot Quaternion
figure;
plot(timeNumeric(1:k), QHistory(:, 1:k));
xlabel('Time [s]');
ylabel('Quaternion Components');
legend('q0', 'q1', 'q2', 'q3');
title('CubeSat Quaternion with Modified B-dot Control');

% Plot Control Torque
figure;
plot(timeNumeric(1:k), controlTorqueHistory(:, 1:k));
xlabel('Time [s]');
ylabel('Control Torque [Nm]');
legend('Tx', 'Ty', 'Tz');
title('Control Torque History');

fprintf('Simulation completed.\n');

%% Function to convert Keplerian elements to ECI coordinates
function [r, v] = keplerian_to_eci(a, ecc, inc, raan, aop, nu, mu)
    % Convert angles to radians
    inc = deg2rad(inc);
    raan = deg2rad(raan);
    aop = deg2rad(aop);
    nu = deg2rad(nu);
    
    % Calculate the distance and speed in the orbital plane
    p = a * (1 - ecc^2);
    r_orbit = p / (1 + ecc * cos(nu));
    v_orbit = sqrt(mu / p);
    
    % Position in the orbital plane
    r_pqw = [r_orbit * cos(nu); r_orbit * sin(nu); 0];
    v_pqw = [-v_orbit * sin(nu); v_orbit * (ecc + cos(nu)); 0];
    
    % Rotation matrices
    R3_W = [cos(raan) -sin(raan) 0; sin(raan) cos(raan) 0; 0 0 1];
        R1_i = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
    R3_w = [cos(aop) -sin(aop) 0; sin(aop) cos(aop) 0; 0 0 1];
    
    % Combined rotation matrix
    Q_pqw_to_eci = R3_W * R1_i * R3_w;
    
    % Position and velocity in ECI frame
    r = Q_pqw_to_eci * r_pqw;
    v = Q_pqw_to_eci * v_pqw;
end


function [lat, lon, alt] = ecef2lla(r)
    % ECEF to latitude, longitude, and altitude conversion
    % Assumes WGS84 ellipsoid
    a = 6378137.0; % Semi-major axis
    f = 1 / 298.257223563; % Flattening
    e2 = 2 * f - f^2; % Square of eccentricity

    x = r(1);
    y = r(2);
    z = r(3);

    lon = atan2(y, x);
    
    % Iterative solution
    p = sqrt(x^2 + y^2);
    lat = atan2(z, p * (1 - e2));
    lat_old = 0;
    while abs(lat - lat_old) > 1e-10
        lat_old = lat;
        N = a / sqrt(1 - e2 * sin(lat)^2);
        alt = p / cos(lat) - N;
        lat = atan2(z, p * (1 - e2 * (N / (N + alt))));
    end

    % Convert to degrees
    lat = rad2deg(lat);
    lon = rad2deg(lon);
    alt = alt / 1000; % Convert to km
end

function [r,v] = propagateOrbit(time, rEpoch, vEpoch)
    % Orbital propagation using Keplerian two-body problem for demonstration
    mu = 3.986e14; %[m^3/s^2]
    dt = time(2) - time(1);
    
    % Magnitude of position and velocity vectors
    r_mag = norm(rEpoch);
    v_mag = norm(vEpoch);
    
    % Specific angular momentum
    h = cross(rEpoch, vEpoch);
    h_mag = norm(h);
    
    % Eccentricity vector
    e = (cross(vEpoch, h) / mu) - (rEpoch / r_mag);
    e_mag = norm(e);
    
    % Semi-major axis
    a = 1 / ((2 / r_mag) - (v_mag^2 / mu));
    
    % Mean motion
    n = sqrt(mu / a^3);
    
    % Mean anomaly
    E = acos(dot(e, rEpoch) / (e_mag * r_mag));
    M = E - e_mag * sin(E);
    
    % New mean anomaly
    M_new = M + n * dt;
    
    % Solve Kepler's equation for E_new
    E_new = M_new;
    for i = 1:1000
        E_new = M_new + e_mag * sin(E_new);
    end
    
    % True anomaly
    nu_new = 2 * atan(sqrt((1 + e_mag) / (1 - e_mag)) * tan(E_new / 2));
    
    % Position and velocity in the orbital plane
    r_orbit = a * (1 - e_mag^2) / (1 + e_mag * cos(nu_new));
    r_pqw = [r_orbit * cos(nu_new); r_orbit * sin(nu_new); 0];
    v_pqw = [-sqrt(mu / (a * (1 - e_mag^2))) * sin(nu_new); sqrt(mu / (a * (1 - e_mag^2))) * (e_mag + cos(nu_new)); 0];
    
    % Rotation matrices
    inc = acos(h(3) / h_mag);
    raan = atan2(h(1), -h(2));
    aop = atan2(e(3), e(1));
    
    R3_W = [cos(raan) -sin(raan) 0; sin(raan) cos(raan) 0; 0 0 1];
    R1_i = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
    R3_w = [cos(aop) -sin(aop) 0; sin(aop) cos(aop) 0; 0 0 1];
    
    % Combined rotation matrix
    Q_pqw_to_eci = R3_W * R1_i * R3_w;
    
    % Position and velocity in ECI frame
    r = Q_pqw_to_eci * r_pqw;
    v = Q_pqw_to_eci * v_pqw;
    
    % Ensure r and v are column vectors
    r = r(:);
    v = v(:);
end
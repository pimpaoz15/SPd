% This script runs just before the Simulink model begins to propagate the
% attitude to define constants. Internal models can be turned off here.
%% -----------------------------------------------------------------------------
%   Copyright (c) 2010-2018 Samir A. Rawashdeh
%   Electrical and Computer Engineering
%   University of Michigan - Dearborn
%  
%   All rights reserved. 
%   
%   Redistribution and use in source and binary forms, with or without 
%   modification, are permitted provided that the following conditions are 
%   met:
%   
%       * Redistributions of source code must retain the above copyright 
%         notice, this list of conditions and the following disclaimer.
%       * Redistributions in binary form must reproduce the above copyright 
%         notice, this list of conditions and the following disclaimer in 
%         the documentation and/or other materials provided with the distribution
%         
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%   POSSIBILITY OF SUCH DAMAGE.
%  
% ----------------------------------------------------------------------------
clc
clear all

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

%% Boira parameters

boira_mass = 2; % Boira's mass (kg)
%boira_Imax = [0.1 0.1 0.15]; %Maximum currents to circulate through the x, y and z magnetorquers (A)

% Maximum magnetic torque that can be produced by the torquers and the coil
boira_max_magnetic_torque = [0.2 0.2 0.12]; % (A m^2) (body frame)
boira_magnetic_gains = [2.3 2.3 2.1]; % Magnetic gains of the manetorquers (A*m²/A)
%boira_magnets = boira_magnetic_gains.*boira_Imax;
boira_inertia = [0.0046 0 0; 0 0.0046 0; 0 0 0.00145]; % Boira's inertia tensor (kg*m²)
boira_V_hyst = [0 0 0]; % Boira's histeresis voltages (unknown) (V)

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

% Initial desired magnetic moment components (must be defined before
% simulation but have no real meaning):
m0 = [0 0 0];

%% PID constants

% Kp_pid = [-10E-11 -10E-11 -10E-11 -10E-11];
% Ki_pid = [-10E-07 -10E-07 -10E-07 -10E-07];
% Kd_pid = [-10E-04 -10E-04 -10E-04 -10E-04];

%%Model with no gravity and no mm models without tuning
% Kp_pid1 = -0.1;
% Ki_pid1 = 0;
% Kd_pid1 = -0.01;
% 
% %%Model with no gravity and no mm models
% Kp_pid1 = -0.05;
% Ki_pid1 = -0.0000005;
% Kd_pid1 = -0.6;
% 
% %%Model with gravity and no mm models
% Kp_pid1 = -10e-5;
% Ki_pid1 = -10e-12;
% Kd_pid1 = -1.1*10e-3;

%%Model with gravity and mm models
Kp_pid1 = -100000;
Ki_pid1 = -100000;
Kd_pid1 = -100000;

Kp_pid = [Kp_pid1 Kp_pid1 Kp_pid1 Kp_pid1];
Ki_pid = [Ki_pid1 Ki_pid1 Ki_pid1 Ki_pid1];
Kd_pid = [Kd_pid1 Kd_pid1 Kd_pid1 Kd_pid1];



%%Model 3
% Kp_pid = -10E-11;
% Ki_pid = -10E-07;
% Kd_pid = -10E-04;

%% Load parameters from GUI
load sim_parameters.mat
SNAP_sim_params = rmfield(SNAP_sim_params, 'name');

%% Enables
EN_2BodyModel = 1;
EN_GG = 1;
EN_Hyst = 0;
EN_Mag = 1;
EN_Aero = SNAP_sim_params.aero_enabled;  %value taken from checkbox in GUI

% display summary to terminal
{'Two-body Model', EN_2BodyModel;
    'Gravity Gradient', EN_GG;
    'Magnetic Hysteresis', EN_Hyst;
    'Magnetic Dipole', EN_Mag;
    'Aerodynamics', EN_Aero}
disp('The above variables indicate which models are enabled and will affect the simulation')


%% Simulation length
SNAP_sim_length = SNAP_sim_params.sim_length;

%% Julian Date
SNAP_JD = SNAP_sim_params.JD;
SNAP_starting_ECI2ECEF_dcm =dcmeci2ecef('IAU-2000/2006',SNAP_sim_params.DateVector)

% Possible override:
%SNAP.JD = mjuliandate(2024,10,10,12,21,0); %Calculate Julian date for October 10, 2024, at 12:21:00 p.m.

%% Orbit Parameters, initial position and velocity
initial_velocity = SNAP_sim_params.initial_velocity ;
initial_position = SNAP_sim_params.initial_position ;

% Possible override (calculation based on altitude and inclination only)
%orbit_alt = 500.e3;
%inclination = 90* pi/180;
%initial_position = [0 (orbit_alt+SNAP_CONST_r_earth) 0];
%initial_velocity = [cos(inclination) 0 sin(inclination)] * 7683;
%initial_velocity = [0 0 0];

%% Mass, inertia, magnet strength, hysteresis material volume

SNAP_sat_mass = boira_mass;
SNAP_sat_inertia = boira_inertia;
% SNAP_magnets = SNAP_sim_params.magnets;
SNAP_V_hyst = boira_V_hyst;

SNAP_Hc = SNAP_sim_params.Hc;
SNAP_Br = SNAP_sim_params.Br;
SNAP_Bs = SNAP_sim_params.Bs;

%% Initial rotation rate
SNAP_init_rates = SNAP_sim_params.init_rates;

%% Constants
SNAP_CONST_m_earth = 5.9736e24;
SNAP_CONST_G_earth = 6.6742e-11;
SNAP_CONST_r_earth = 6.36101e6;
SNAP_CONST_mu0 = 4*pi*1e-7;

%Earth Mag Model
SNAP_CONST_earth_mag_dipole = [sin(168.6*pi/180)* cos(109.3*pi/180), sin(168.6*pi/180)* sin(109.3*pi/180), cos(168.6*pi/180)] ; %unit vector, in ECEF
SNAP_CONST_a3H0 = 7.943e15; % Wb.m


%% Aerodynamic model parameters
if EN_Aero == 1
    load("SNAP_Basic_2U.mat");
    
else  % dummy values so that aero block in simulink model does not complain about missing variables
    SNAP_aeromodel.T = zeros (37,32);
    SNAP_aeromodel.pitch = 1:37;
    SNAP_aeromodel.roll = 1:32;
    SNAP_aeromodel.av_density_vs_alt = zeros(1,11);
    SNAP_aeromodel.alt_range = 1:11;
        
end


SNAPaero.T = SNAP_aeromodel.T;
SNAPaero.pitch = SNAP_aeromodel.pitch;
SNAPaero.roll = SNAP_aeromodel.roll;
SNAPaero.av_density_vs_alt = SNAP_aeromodel.av_density_vs_alt;
SNAPaero.alt_range = SNAP_aeromodel.alt_range;


%% 

% % %%Plots
% figure()
% subplot(3,1,1);
% plot(out_angs2nadir.signals.values(:,1))
% title("x-axis")
% subplot(3,1,2);
% plot(out_angs2nadir.signals.values(:,2))
% title("y-axis")
% subplot(3,1,3);
% plot(out_angs2nadir.signals.values(:,3))
% title("z-axis")
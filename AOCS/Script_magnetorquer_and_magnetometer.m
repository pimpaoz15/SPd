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


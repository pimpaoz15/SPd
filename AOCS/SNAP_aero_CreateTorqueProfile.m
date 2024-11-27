%% This script generates torque profiles for SNAP's aerodynamic model.
%
% -----------------------------------------------------------------------------
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

clear
addpath('libaero')


%% Define "home" volume
Resolution = 0.5; % cm cube per dot
res = 360;  %home volume dimentions
home_volume = zeros(res,res,res) * NaN;

%% Define satellite Shape
[origin, home_volume] = draw_cubesat(30, 0, home_volume, Resolution);
origin  = [res res res]/2;  % overwrite, numerical rounding caused odd results.

SNAP_aeromodel.PointCloudModel = home_volume;

plot_volume(home_volume)
drawnow

%% rotate volume
roll = (0:3:45) * pi/180;
pitch = (0:5:180) * pi/180;
yaw = 0 * pi/180;

T = zeros(length(roll), length(pitch));
for iroll = 1:length(roll)
    for ipitch = 1:length(pitch)
        
        DCM_roll = angle2dcm( roll(iroll), 0 , 0, 'XYZ');
        DCM_pitch = angle2dcm( 0, pitch(ipitch) , 0, 'XYZ');
        DCM = DCM_roll*DCM_pitch;
        rot_volume = rotate_volume(home_volume, DCM, origin);
        
        % plot_volume(rot_volume)
        
        % Torque in body frame
        temp = calc_torque(rot_volume, origin);
        
        T(iroll, ipitch) =  sqrt(temp(1)^2+temp(2)^2+temp(3)^2);
        disp(['pitch: ' num2str(pitch(ipitch)*180/pi) '/180, roll: '  num2str(roll(iroll)*180/pi) '/45' ', Torque =' num2str(temp)])
        
        %           pause
    end
    
    %         plot(T)
    
end
%T(T<1e-8) = 0;

SNAP_aeromodel.T = [T' fliplr(T')];

SNAP_aeromodel.pitch = pitch;
SNAP_aeromodel.roll = [roll roll+roll(length(roll))];

SNAP_aeromodel.alt_range = (200:50:700); % km, altitudes
SNAP_aeromodel.lo_density_vs_alt = [1.78e-10 3.35e-11 8.19e-12 2.34e-12 7.32e-13 2.47e-13 8.98e-14 3.63e-14 1.68e-14 9.14e-15 5.74e-15]; %% averages, Kg/m3
SNAP_aeromodel.av_density_vs_alt = [2.53e-10 6.24e-11 1.95e-11 6.98e-12 2.72e-12 1.13e-12 4.89e-13 2.21e-13 1.04e-13 5.15e-14 2.71e-14]; %% averages, Kg/m3
SNAP_aeromodel.up_density_vs_alt = [3.52e-10 1.06e-10 3.96e-11 1.66e-11 7.55e-12 3.61e-12 1.8e-12 9.25e-13 4.89e-13 2.64e-13 1.47e-13 ]; %% averages, Kg/m3

save('SNAP_Basic_3U','SNAP_aeromodel')




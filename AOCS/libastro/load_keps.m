% [R V Epoch JD Name] = load_keps(file)
%
% Translates TLEs to orbit initial conditions (position, velocity and
% epoch). If the file contains multiple TLEs, the first one is used.
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


function [R V Epoch Name] = load_keps(file)

sat = read_tloes(file);
i=1;
% for i = 1:length(sat)
    
%% calculating semi major axis
n = sat(i).mean_motion*2*pi/(24*60*60);
mu = 5.9736e24 * 6.6742e-11;
a = (mu/n^2)^(1/3); 

%% calculating true anomaly
nu = nuFromM(sat(i).mean_anomaly*pi/180,sat(i).eccentricity);

%% calculate R and V in km/sec
[R,V] = randv(a/1000, ...
    sat(i).eccentricity, ...
    sat(i).inclination*pi/180, ...
    sat(i).RAAN*pi/180, ...
    sat(i).argp*pi/180, ...
    nu);

% in m/s
R = R*1000;
V = V*1000;

Epoch = sat(i).epoch.year + sat(1).epoch.day/365.2425;

Name = sat(i).name;

%  disp([sat(i).name ', alt:' num2str((a-6357000)/1000) ', inc:' num2str(sat(i).inclination)])
%  alt = (sqrt(R(1)^2+R(2)^2+R(3)^2)-6357000)/1000
%  vel = (sqrt(V(1)^2+V(2)^2+V(3)^2))/1000
%  
 % end

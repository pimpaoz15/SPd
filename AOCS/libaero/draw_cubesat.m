% Function to create the point cloud of a basic CubeSat. 
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

function [cg, volume] = draw_cubesat(body_length, cg_offset, volume, Resolution)

volume = volume * NaN;
len = body_length / Resolution;  % length in dots = 30cm * Resolution
wid = 10 / Resolution;  % width in dots = 10cm * Resolution
cg = round(size(volume)/2); % temporary cg - center of the volume



%% main satellite body
x_range = (cg-len/2):(cg+len/2);
y_range = (cg-wid/2):(cg+wid/2);
z_range = (cg-wid/2):(cg+wid/2);
volume(x_range,y_range,z_range) = 1;

x_range = (cg-len/2)+2:(cg+len/2)-2;
y_range = (cg-wid/2)+2:(cg+wid/2)-2;
z_range = (cg-wid/2)+2:(cg+wid/2)-2;
volume(x_range,y_range,z_range) = NaN;

%% Example, to add deployables, like antennae or booms:
% x_range = 200:220;
% y_range = 189:190;
% z_range = 189;
% volume(x_range,y_range,z_range) = 1;



%% CG Correction
%% IMPORTANT:  Best to hardcode CG outside of this function, variable called "origin"
[x y z] = find_1s_in_volume(volume);
cg = round([sum(x)/length(x), sum(y)/length(y), sum(z)/length(z)]);  % assumes satellite is uniformly distributed
cg(1) = cg(1) + cg_offset/Resolution;  


        
       




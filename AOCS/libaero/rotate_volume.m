%% rotates the volume containing the satelltie diefined as '1s', by the
% direction cosine matrix (DCM) about the 'origin'
% takes 0.8 seconds
%
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

function rotated_volume = rotate_volume(volume, DCM, origin)

rotated_volume = NaN(size(volume)); % initialize rotated volume to the right size


%% find x y z of satellite surfaces
[x y z] = find_1s_in_volume(volume);

%% rotate each volume element by the DCM around the origin, note the new
%% indecies, and set them to 1 in the new volume.

for i = 1:length(x)  % loop through all points in volume
    rot =  ([x(i),y(i),z(i)] - origin) * DCM;
    rot = rot + origin;
%     rot(rot<1) = 1;
%     rot(rot>length(volume)) = length(volume);
    
    rot = round(rot);
    rotated_volume(rot(1),rot(2),rot(3)) = 1;
    

    
end



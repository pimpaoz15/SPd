% Function to plot a point cloud for debugging purposes.
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

function plot_volume(volume)

% figure
% hold on

[x y z] = find_1s_in_volume(volume);
plot3(x,y,z,'.k')
xlabel('x')
ylabel('y')
zlabel('z')

grid on
axis([1 size(volume,1) 1 size(volume,2) 1 size(volume,3)])

% 
% [X,Y,Z] = meshgrid(1:size(volume,1), 1:size(volume,2), 1:size(volume,3));
% 
% SLICE(X,Y,Z,volume,x,y,z)





% function plot_volume(volume)
% 
% figure
% hold on
% 
% [x y z] = find_1s_in_volume(volume);
% for i = x(1):x(end)
%     for j = y(1):y(end)
%         for k = z(1):z(end)
% 
% plot3(i,j,k,'.','Color',[i/80,j/80,k/80])
% xlabel('x')
% ylabel('y')
% zlabel('z')
%         end
%     end
% end
% 
% grid on
% axis([1 size(volume,1) 1 size(volume,2) 1 size(volume,3)])
% 
% % 
% % [X,Y,Z] = meshgrid(1:size(volume,1), 1:size(volume,2), 1:size(volume,3));
% % 
% % SLICE(X,Y,Z,volume,x,y,z)
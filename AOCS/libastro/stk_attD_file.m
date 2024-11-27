% stk_attD_file(ATT,epoch,name)
%
% Generates attitude file for STK, given the variable "ATT" containing the 
% DCM and time vectors as a 'simout' from Simulink. 
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

function stk_attD_file(ATT,epoch,name)


file = fopen(name,'w');

fprintf(file,'stk.v.5.0\n');
fprintf(file,'BEGIN Attitude\n\n');

fprintf(file,['NumberOfAttitudePoints  ' num2str(length(ATT.signals.values)) ' \n']);
fprintf(file,'BlockingFactor          20 \n');
fprintf(file,'InterpolationOrder      1 \n');
fprintf(file,'CentralBody             Earth \n');
fprintf(file,'CoordinateAxes        	J2000 \n\n');

fprintf(file,'AttitudeTimeDCM \n\n');

DCM = ATT.signals.values(:,:,:);

for i = 1:length(ATT.signals.values)
    fprintf(file,'%8.12E %f %f %f %f %f %f %f %f %f\n', ...
        ATT.time(i), DCM(1,1,i), DCM(1,2,i), DCM(1,3,i),...
        DCM(2,1,i), DCM(2,2,i), DCM(2,3,i), ...
        DCM(3,1,i), DCM(3,2,i), DCM(3,3,i)  );
end


fprintf(file,'\nEND Attitude\n\n');
fclose(file);

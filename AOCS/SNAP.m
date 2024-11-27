% This is the main script for the Smart Nano-satellite Attitude Propagator.
% Running this m-file will bring up the GUI
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

function varargout = SNAP(varargin)

%
% Edit the above text to modify the response to help SNAP
% Last Modified by GUIDE v2.5 15-Aug-2018 14:14:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SNAP_OpeningFcn, ...
    'gui_OutputFcn',  @SNAP_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

end

%% --- Executes just before SNAP is made visible.
function SNAP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SNAP (see VARARGIN)

addpath('libastro')

% Choose default command line output for SNAP
handles.output = hObject;

handles.status.sim = 0;

% Update handles structure
guidata(hObject, handles);

% Plot Map
axes(handles.axes1)
earth = imread('earth.jpg');
lv= size(earth,1);
lh= size(earth,2);

lats =  (1:lv)*180/lv - 90;
lons =  (1:lh)*360/lh - 180;

image(lons, -lats, earth(1024:-1:1,:,:))
grid on
set(handles.axes1,'XTick',[-180:30:180]);
set(handles.axes1,'YTick',[-90:30:90]);
set(handles.axes1,'YTickLabel',[90:-30:-90]);
set(handles.axes1,'Xcolor',0.3*ones(1,3));
set(handles.axes1,'Ycolor',0.3*ones(1,3));
title('Ground Track')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')


%% Load default parameters
load sim_parameters
writeGUIfields(SNAP_sim_params, hObject, eventdata, handles)


% UIWAIT makes SNAP wait for user response (see UIRESUME)
% uiwait(handles.figure1);

end


%% POP MENU --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)

item = get(hObject,'Value') %returns selected item from popupmenu1
hold off

if (get(handles.checkbox4,'Value') == 1)
    figure
else
    axes(handles.axes1); %set the current axes to axes1
end

switch(item)
    case 1 %HEADING
        
    case 2 %point cloud
        if (get(handles.checkbox5,'Value') == 0)
            msgbox('Aero model checkbox is disabled.')
        else
            load(get(handles.edit34,'string'));
        end
        
        
        if ~exist('SNAP_aeromodel','var')
            msgbox('Error: expected variables not found.')
        else
            plot_volume(SNAP_aeromodel.PointCloudModel)
            view(-50,40)
            zoom(5.5)
        end
        
    case 3 %1D aero torque profile
        if (get(handles.checkbox5,'Value') == 0)
            msgbox('Aero model checkbox is disabled.')
        else
            load(get(handles.edit34,'string'));
        end
        
        
        if ~exist('SNAP_aeromodel','var')
            msgbox('Error: expected variables not found.')
        else
            y = [-SNAP_aeromodel.T(:,16); fliplr(SNAP_aeromodel.T(:,16))];
            x = [fliplr(SNAP_aeromodel.pitch) -SNAP_aeromodel.pitch]*180/pi;
            plot(x,y)
            set(gca,'XTick',[-180:30:180])
            set(gca,'XLim',[-180,180])
            grid on
            xlabel('Pitch Angles (degrees)')
            ylabel('Pitch Torque - Normalized (N.m / [Velocity^2 . Atmospheric Density)]')
            
        end
        
        
    case 4 %2D mesh plot of aero torque profile
        if (get(handles.checkbox5,'Value') == 0)
            msgbox('Aero model checkbox is disabled.')
        else
            load(get(handles.edit34,'string'));
        end
        
        
        if ~exist('SNAP_aeromodel','var')
            msgbox('Error: expected variables not found.')
        else
            
            mesh(SNAP_aeromodel.roll*180/pi, ...
                SNAP_aeromodel.pitch*180/pi,...
                SNAP_aeromodel.T)
            set(gca,'XTick',[0:30:90]);
            set(gca,'YTick',[0:30:180]);
            title('Aerodynamic torque profile')
            ylabel('Pitch Angle (degrees)')
            xlabel('Roll Angle (degrees)')
            zlabel('Pitch Torque Factor (N.m / [Velocity^2 * Air Density)')
            
        end
        
    case 5  %HEADING
        
        
    case 6 %ground trace
        if(handles.status.sim == 0)
            msgbox('No data to plot, run a simulation first.')
        else
            
            earth = imread('earth.jpg');
            lv= size(earth,1);
            lh= size(earth,2);
            
            lats =  (1:lv)*180/lv - 90;
            lons =  (1:lh)*360/lh - 180;
            
            image(lons, -lats, earth(1024:-1:1,:,:))
            grid on
            set(handles.axes1,'XTick',[-180:30:180]);
            set(handles.axes1,'YTick',[-90:30:90]);
            set(handles.axes1,'YTickLabel',[90:-30:-90]);
            set(handles.axes1,'Xcolor',0.3*ones(1,3));
            set(handles.axes1,'Ycolor',0.3*ones(1,3));
            title('Ground Track')
            xlabel('Longitude (degrees)')
            ylabel('Latitude (degrees)')
            
            hold on
            lla = ecef2lla2(handles.out.position);
            x_lla = lla(1:100:end,2);
            y_lla = -lla(1:100:end,1);
            
            plot(x_lla,y_lla,'.y');
        end
        
    case 7 %euler angles
        if(handles.status.sim == 0)
            msgbox('No data to plot, run a simulation first.')
        else
            plot(handles.out.time/60,handles.out.euler*180/pi)
            grid on
            title('Euler Angles')
            xlabel('Time (minutes)')
            ylabel('Angles (degrees)')
            legend('About Roll Axis','About Pitch Axis','About Yaw Axis')
        end
        
    case 8 %angs2mf
        if(handles.status.sim == 0)
            msgbox('No data to plot, run a simulation first.')
        else
            plot(handles.out.time/60,handles.out.angs2mf)
            grid on
            title('Satellite Main Axes Relative to Magnetic Field Vector')
            xlabel('Time (minutes)')
            ylabel('Angles (degrees)')
            legend('Roll Axis','Pitch Axis','Yaw Axis')
        end
        
    case 9 %angs2gg
        if(handles.status.sim == 0)
            msgbox('No data to plot, run a simulation first.')
        else
            plot(handles.out.time/60,handles.out.angs2nadir)
            grid on
            title('Satellite Main Axes Relative to Nadir Vector')
            xlabel('Time (minutes)')
            ylabel('Angles (degrees)')
            legend('Roll Axis','Pitch Axis','Yaw Axis')
        end
        
    case 10 %angs2vv
        if(handles.status.sim == 0)
            msgbox('No data to plot, run a simulation first.')
        else
            plot(handles.out.time/60,handles.out.angs2velocity)
            grid on
            title('Satellite Main Axes Relative to Velocity Vector')
            xlabel('Time (minutes)')
            ylabel('Angles (degrees)')
            legend('Roll Axis','Pitch Axis','Yaw Axis')
        end
        
    case 11 %rotation rates
        if(handles.status.sim == 0)
            msgbox('No data to plot, run a simulation first.')
        else
            plot(handles.out.time/60,handles.out.w*180/pi)
            grid on
            title('Satellite Angular Rotation Rates')
            xlabel('Time (minutes)')
            ylabel('Rotation Rates (degrees/sec)')
            legend('About Roll Axis','About Pitch Axis','About Yaw Axis')
        end
        
    case 12 %gg_trq
        if(handles.status.sim == 0)
            msgbox('No data to plot, run a simulation first.')
        else
            plot(handles.out.time/60,handles.out.trq_gg)
            grid on
            title('Gravity Gradient Torque')
            xlabel('Time (minutes)')
            ylabel('Torque (N.m)')
            legend('About Roll Axis','About Pitch Axis','About Yaw Axis')
        end
        
    case 13 %mag_trq
        if(handles.status.sim == 0)
            msgbox('No data to plot, run a simulation first.')
        else
            plot(handles.out.time/60,handles.out.trq_mag)
            grid on
            title('Torque due to Permanent Magnet')
            xlabel('Time (minutes)')
            ylabel('Torque (N.m)')
            legend('About Roll Axis','About Pitch Axis','About Yaw Axis')
        end
        
    case 14 %Hyst_trq
        if(handles.status.sim == 0)
            msgbox('No data to plot, run a simulation first.')
        else
            plot(handles.out.time/60,handles.out.trq_hyst)
            grid on
            title('Torque due to Hysteresis Material')
            xlabel('Time (minutes)')
            ylabel('Torque (N.m)')
            legend('About Roll Axis','About Pitch Axis','About Yaw Axis')
        end
        
    case 15 %Aero_trq
        if(handles.status.sim == 0)
            msgbox('No data to plot, run a simulation first.')
        else
            plot(handles.out.time/60,handles.out.trq_aero)
            grid on
            title('Torque due to Aerodynamics')
            xlabel('Time (minutes)')
            ylabel('Torque (N.m)')
            legend('About Roll Axis','About Pitch Axis','About Yaw Axis')
        end
        
    case 16 %All environmental torques
        if(handles.status.sim == 0)
            msgbox('No data to plot, run a simulation first.')
        else
            
            for i=1:length(handles.out.time)
                tq_aero_abs(i) = sqrt( handles.out.trq_aero(1,i)^2 + handles.out.trq_aero(2,i)^2 + handles.out.trq_aero(3,i)^2 );
                tq_gg_abs(i) = sqrt( handles.out.trq_gg(1,i)^2 + handles.out.trq_gg(2,i)^2 + handles.out.trq_gg(3,i)^2 );
                tq_hyst_abs(i) = sqrt( handles.out.trq_hyst(1,i)^2 + handles.out.trq_hyst(2,i)^2 + handles.out.trq_hyst(3,i)^2 );
                tq_mag_abs(i) = sqrt( handles.out.trq_mag(1,i)^2 + handles.out.trq_mag(2,i)^2 + handles.out.trq_mag(3,i)^2 );
                
            end
            
            time = handles.out.time;
            
            semilogy(time,tq_mag_abs,'b',time,tq_aero_abs,'g--',time,tq_gg_abs,'r.',time,tq_hyst_abs,'m.-');
            %plot(time,tq_mag_abs,'b',time,tq_aero_abs,'g--',time,tq_gg_abs,'r.',time,tq_hyst_abs,'m.-');
            
            title('Magnitude of environmental torques')
            xlabel('Time (hours)')
            ylabel('Torque (N/m^2)')
            legend('Magnetic','Aerodynamic','Gravity Gradient','Magnetic Hysteresis')
            
        end
end

end


%% PUSH BUTTON: Run Simulation
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Read Parameters
SNAP_sim_params = readGUIfields(hObject, eventdata, handles);

%% Pass variables to Simulink Model, and start simulation
save sim_parameters SNAP_sim_params

%% Plot map in preparation for ground track after simulation is done. It's a good indicator that the simulation has completed

% Plot Map
axes(handles.axes1)
earth = imread('earth.jpg');
lv= size(earth,1);
lh= size(earth,2);

lats =  (1:lv)*180/lv - 90;
lons =  (1:lh)*360/lh - 180;

image(lons, -lats, earth(1024:-1:1,:,:))
grid on
set(handles.axes1,'XTick',[-180:30:180]);
set(handles.axes1,'YTick',[-90:30:90]);
set(handles.axes1,'YTickLabel',[90:-30:-90]);
set(handles.axes1,'Xcolor',0.3*ones(1,3));
set(handles.axes1,'Ycolor',0.3*ones(1,3));
title('Ground Track')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
drawnow

%% Run the simulation

% sim('SNAP_Simulink_model_v10_NoAero.slx')
sim('SNAP_Simulink_model_v10.slx')
msgbox('Simulation completed, plotting ground track next.')

%% Sim completed;  extract output variables
handles.out.angs2mf = out_angs2mf.signals.values;
handles.out.angs2velocity = out_angs2velocity.signals.values;
handles.out.angs2nadir = out_angs2nadir.signals.values;
handles.out.euler = out_euler.signals.values;
handles.out.position = out_position.signals.values(:,:);
handles.out.velocity = out_velocity.signals.values(:,:);
handles.out.trq_gg = out_trq_gg.signals.values(:,:);
handles.out.trq_hyst = out_trq_hyst.signals.values(:,:);
handles.out.trq_mag = out_trq_mag.signals.values(:,:);
handles.out.trq_aero = out_trq_aero.signals.values(:,:);
handles.out.w = out_w.signals.values;
handles.out.time = tout;

% signal that a simulation has been run
handles.status.sim = 1;


% plot ground track
axes(handles.axes1)
earth = imread('earth.jpg');
lv= size(earth,1); lh= size(earth,2);
lats =  (1:lv)*180/lv - 90; lons =  (1:lh)*360/lh - 180;
image(lons, -lats, earth(1024:-1:1,:,:))
grid on
set(handles.axes1,'XTick',[-180:30:180]);
set(handles.axes1,'YTick',[-90:30:90]);
set(handles.axes1,'YTickLabel',[90:-30:-90]);
set(handles.axes1,'Xcolor',0.3*ones(1,3));
set(handles.axes1,'Ycolor',0.3*ones(1,3));
title('Ground Track')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')

hold on
lla = ecef2lla2(handles.out.position);
x_lla = lla(1:100:end,2);
y_lla = -lla(1:100:end,1);

plot(x_lla,y_lla,'.y');


% if box is checked, save MAT file
if (get(handles.checkbox2,'Value') == 1)
    sim_results = handles.out;
    sim_results.SNAP_sim_params = SNAP_sim_params;
    filename = ['OUTPUT_' SNAP_sim_params.name];
    save(filename, 'sim_results')
end

% if box is checked, generate STK Attitude File
if (get(handles.checkbox1,'Value') == 1)
    filename = ['OUTPUT_STK_' SNAP_sim_params.name '.a'];
    stk_attD_file(Attitude_DCM, SNAP_sim_params.epoch, filename);
end


% updated handles to make variables global
guidata(hObject, handles);

end


%% Load Orbital info from TLEs
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uigetfile('*.txt','Select the TLEs text file');
if FileName ~= 0
    SNAP_sim_params = readGUIfields(hObject, eventdata, handles);
    [SNAP_sim_params.initial_position SNAP_sim_params.initial_velocity ...
        SNAP_sim_params.epoch SNAP_sim_params.name] = load_keps([PathName FileName]);
    writeGUIfields(SNAP_sim_params, hObject, eventdata, handles);
end

end

%% Write GUI Fields
function writeGUIfields(SNAP_sim_params, hObject, eventdata, handles)
set(handles.edit1,'string',num2str(SNAP_sim_params.initial_position(1)/1000))
set(handles.edit20,'string',num2str(SNAP_sim_params.initial_position(2)/1000))
set(handles.edit21,'string',num2str(SNAP_sim_params.initial_position(3)/1000))

set(handles.edit23,'string',num2str(SNAP_sim_params.initial_velocity(1)/1000))
set(handles.edit24,'string',num2str(SNAP_sim_params.initial_velocity(2)/1000))
set(handles.edit25,'string',num2str(SNAP_sim_params.initial_velocity(3)/1000))

set(handles.edit22,'string',num2str(SNAP_sim_params.m_satellite))

set(handles.edit4,'string',num2str(SNAP_sim_params.inertia_satellite(1,1)))
set(handles.edit6,'string',num2str(SNAP_sim_params.inertia_satellite(1,2)))
set(handles.edit7,'string',num2str(SNAP_sim_params.inertia_satellite(1,3)))
set(handles.edit8,'string',num2str(SNAP_sim_params.inertia_satellite(2,1)))
set(handles.edit9,'string',num2str(SNAP_sim_params.inertia_satellite(2,2)))
set(handles.edit10,'string',num2str(SNAP_sim_params.inertia_satellite(2,3)))
set(handles.edit11,'string',num2str(SNAP_sim_params.inertia_satellite(3,1)))
set(handles.edit12,'string',num2str(SNAP_sim_params.inertia_satellite(3,2)))
set(handles.edit13,'string',num2str(SNAP_sim_params.inertia_satellite(3,3)))

set(handles.edit14,'string',num2str(SNAP_sim_params.magnets(1)))
set(handles.edit15,'string',num2str(SNAP_sim_params.magnets(2)))
set(handles.edit16,'string',num2str(SNAP_sim_params.magnets(3)))

set(handles.edit17,'string',num2str(SNAP_sim_params.V_hyst(1)*1e6))
set(handles.edit18,'string',num2str(SNAP_sim_params.V_hyst(2)*1e6))
set(handles.edit19,'string',num2str(SNAP_sim_params.V_hyst(3)*1e6))

set(handles.edit26,'string',num2str(SNAP_sim_params.Hc))
set(handles.edit27,'string',num2str(SNAP_sim_params.Br))
set(handles.edit28,'string',num2str(SNAP_sim_params.Bs))

set(handles.edit29,'string', num2str(SNAP_sim_params.sim_length/60/60));

set(handles.edit30,'string', num2str(180/pi*SNAP_sim_params.init_rates(1)));
set(handles.edit31,'string', num2str(180/pi*SNAP_sim_params.init_rates(2)));
set(handles.edit32,'string', num2str(180/pi*SNAP_sim_params.init_rates(3)));

set(handles.edit3,'string', num2str(SNAP_sim_params.epoch));

set(handles.edit33,'string',SNAP_sim_params.name);

set(handles.edit34,'string',SNAP_sim_params.aero_file);
set(handles.checkbox5,'Value', SNAP_sim_params.aero_enabled);

end


%% Read GUI Fields
function SNAP_sim_params = readGUIfields(hObject, eventdata, handles)
% [SNAP_sim_params.initial_position, SNAP_sim_params.initial_velocity] = load_keps('TLEs_DelfiC3.txt');

SNAP_sim_params.initial_position(1,1) = 1000*str2num(get(handles.edit1,'string'));
SNAP_sim_params.initial_position(1,2) = 1000*str2num(get(handles.edit20,'string'));
SNAP_sim_params.initial_position(1,3) = 1000*str2num(get(handles.edit21,'string'));

SNAP_sim_params.initial_velocity(1,1) = 1000*str2num(get(handles.edit23,'string'));
SNAP_sim_params.initial_velocity(1,2) = 1000*str2num(get(handles.edit24,'string'));
SNAP_sim_params.initial_velocity(1,3) = 1000*str2num(get(handles.edit25,'string'));

SNAP_sim_params.m_satellite = str2num(get(handles.edit22,'string'));

SNAP_sim_params.inertia_satellite(1,1) = str2num(get(handles.edit4,'string'));
SNAP_sim_params.inertia_satellite(1,2) = str2num(get(handles.edit6,'string'));
SNAP_sim_params.inertia_satellite(1,3) = str2num(get(handles.edit7,'string'));
SNAP_sim_params.inertia_satellite(2,1) = str2num(get(handles.edit8,'string'));
SNAP_sim_params.inertia_satellite(2,2) = str2num(get(handles.edit9,'string'));
SNAP_sim_params.inertia_satellite(2,3) = str2num(get(handles.edit10,'string'));
SNAP_sim_params.inertia_satellite(3,1) = str2num(get(handles.edit11,'string'));
SNAP_sim_params.inertia_satellite(3,2) = str2num(get(handles.edit12,'string'));
SNAP_sim_params.inertia_satellite(3,3) = str2num(get(handles.edit13,'string'));

SNAP_sim_params.magnets(1) = str2num(get(handles.edit14,'string'));
SNAP_sim_params.magnets(2) = str2num(get(handles.edit15,'string'));
SNAP_sim_params.magnets(3) = str2num(get(handles.edit16,'string'));

SNAP_sim_params.V_hyst(1)  = 1e-6* str2num(get(handles.edit17,'string'));
SNAP_sim_params.V_hyst(2)  = 1e-6* str2num(get(handles.edit18,'string'));
SNAP_sim_params.V_hyst(3)  = 1e-6* str2num(get(handles.edit19,'string'));

SNAP_sim_params.Hc = str2num(get(handles.edit26,'string'));
SNAP_sim_params.Br = str2num(get(handles.edit27,'string'));
SNAP_sim_params.Bs = str2num(get(handles.edit28,'string'));

SNAP_sim_params.sim_length = 60 * 60 * eval(get(handles.edit29,'string'));

SNAP_sim_params.init_rates(1) = pi/180 * str2num(get(handles.edit30,'string'));
SNAP_sim_params.init_rates(2) = pi/180 * str2num(get(handles.edit31,'string'));
SNAP_sim_params.init_rates(3) = pi/180 * str2num(get(handles.edit32,'string'));

SNAP_sim_params.epoch = str2num(get(handles.edit3,'string'));
[yy mm dd HH MM SS] = datevec(datenum(fix(SNAP_sim_params.epoch) ...
    ,0,(SNAP_sim_params.epoch - fix(SNAP_sim_params.epoch)) * 365.2425));
SNAP_sim_params.DateVector = [yy mm dd HH MM SS]; 
SNAP_sim_params.JD = mjuliandate([yy mm dd HH MM SS]);


SNAP_sim_params.name = get(handles.edit33,'string');

SNAP_sim_params.aero_file = get(handles.edit34,'string');
SNAP_sim_params.aero_enabled = get(handles.checkbox5,'Value')

end



%% HELP BUTTON
function pushbutton5_Callback(hObject, eventdata, handles)

web help.txt
end


%% UNUSED CALLBACKS AND CREATE FUNCTIONS
function edit29_Callback(hObject, eventdata, handles) end
function edit29_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit30_Callback(hObject, eventdata, handles) end
function edit30_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit31_Callback(hObject, eventdata, handles) end
function edit31_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit32_Callback(hObject, eventdata, handles) end
function edit32_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit1_Callback(hObject, eventdata, handles) end

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit33_Callback(hObject, eventdata, handles) end

function edit33_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit2_Callback(hObject, eventdata, handles) end

function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit3_Callback(hObject, eventdata, handles) end
function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit4_Callback(hObject, eventdata, handles) end
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit6_Callback(hObject, eventdata, handles) end
function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit7_Callback(hObject, eventdata, handles) end
function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit8_Callback(hObject, eventdata, handles) end
function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit9_Callback(hObject, eventdata, handles) end
function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit10_Callback(hObject, eventdata, handles) end
function edit10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit11_Callback(hObject, eventdata, handles) end
function edit11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit12_Callback(hObject, eventdata, handles) end
function edit12_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit13_Callback(hObject, eventdata, handles) end
function edit13_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit14_Callback(hObject, eventdata, handles) end
function edit14_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit15_Callback(hObject, eventdata, handles) end
function edit15_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit16_Callback(hObject, eventdata, handles) end
function edit16_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit17_Callback(hObject, eventdata, handles) end
function edit17_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit18_Callback(hObject, eventdata, handles) end
function edit18_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit19_Callback(hObject, eventdata, handles) end
function edit19_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit22_Callback(hObject, eventdata, handles) end
function edit22_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit26_Callback(hObject, eventdata, handles) end
function edit26_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit27_Callback(hObject, eventdata, handles) end
function edit27_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit28_Callback(hObject, eventdata, handles) end
function edit28_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit20_Callback(hObject, eventdata, handles) end
function edit20_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit21_Callback(hObject, eventdata, handles) end
function edit21_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit23_Callback(hObject, eventdata, handles) end
function edit23_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit24_Callback(hObject, eventdata, handles) end
function edit24_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit25_Callback(hObject, eventdata, handles) end
function edit25_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles) end
function checkbox2_Callback(hObject, eventdata, handles) end
function checkbox4_Callback(hObject, eventdata, handles) end
function checkbox5_Callback(hObject, eventdata, handles) end

function edit34_Callback(hObject, eventdata, handles) end

% --- Outputs from this function are returned to the command line.
function varargout = SNAP_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
end


function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit SNAP_Simulink_Init_Function.m
end

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit SNAP_aero_init.m
end

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (get(handles.checkbox5,'Value') == 0)
    msgbox('Aero model checkbox is disabled. Doing nothing.')
    return
end

[FileName,PathName] = uigetfile('*.mat','Select aerodynamic model variables');
if FileName ~= 0
    SNAP_sim_params = readGUIfields(hObject, eventdata, handles);
    SNAP_sim_params.aero_file = FileName;
    writeGUIfields(SNAP_sim_params, hObject, eventdata, handles);
end

load(SNAP_sim_params.aero_file);

if ~exist('SNAP_aeromodel','var')
    msgbox('Error: loaded aerodynamic model file did not contain expected variables.')
    return
end

y = [-SNAP_aeromodel.T(:,16); fliplr(SNAP_aeromodel.T(:,16))];
x = [fliplr(SNAP_aeromodel.pitch) -SNAP_aeromodel.pitch]*180/pi;
plot(x,y)
set(gca,'XTick',[-180:30:180])
set(gca,'XLim',[-180,180])
grid on
xlabel('Pitch Angles (degrees)')
ylabel('Pitch Torque - Normalized (N.m / [Velocity^2 . Atmospheric Density)]')

end



%% Alternate plot that shows a 2D mesh to include roll angle axis
% mesh(roll*180/pi,pitch*180/pi,T)
% set(gca,'XTick',[0:30:90]);
% set(gca,'YTick',[0:30:180]);
% title('Aerodynamic torque profile')
% ylabel('Pitch Angle (degrees)')
% xlabel('Roll Angle (degrees)')
% zlabel('Pitch Torque Factor (N.m / [Velocity^2 * Air Density)')


% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


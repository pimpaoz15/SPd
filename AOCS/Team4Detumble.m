% Existing MATLAB script with modifications for appending to Word document

% Plot Configuration
reset(groot)
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 21);
set(groot, 'defaultAxesGridAlpha', 0.3);
set(groot, 'defaultAxesLineWidth', 2.25);
set(groot, 'defaultAxesXMinorTick', 'on');
set(groot, 'defaultAxesYMinorTick', 'on');
set(groot, 'defaultFigureRenderer', 'painters');
set(groot, 'defaultLegendBox', 'on');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLegendLocation', 'best');
set(groot, 'defaultLineLineWidth', 2.0);
set(groot, 'defaultLineMarkerSize', 20);
set(groot, 'defaultTextInterpreter', 'latex');

% Initial parameters
a = 396 + 6378;
mu = 398600;
nu = 2*pi*(a^3/mu)^0.5;

W_val = out_w.signals.values;
% Bdotvals = out_bdot.signals.values;
% 
% [~, ~, dim3] = size(Bdotvals);
% reshaped_array = reshape(Bdotvals, 3, dim3);
% Bdotx = reshaped_array(1,:);
% Bdoty = reshaped_array(2,:);
% Bdotz = reshaped_array(3,:);

wx = W_val(:,1)*180/pi;
wy = W_val(:,2)*180/pi;
wz = W_val(:,3)*180/pi;

t = out_w.time;
OrbNum = t/nu;

% Plot Angular rate dampening
figure;
plot(t, wx)
hold on
plot(t, wy)
hold on
plot(t,wz)
title("Angular rate dampening")
xlabel("Time (s)")
ylabel("Angular rate (rad/s)")
legend("w_x","w_y","w_z")
grid on
saveas(gcf, "Images/Detumble_s.png")

figure;
plot(OrbNum, wx)
hold on
plot(OrbNum, wy)
hold on
plot(OrbNum,wz)
title("Angular rate dampening")
xlabel("Orbit period number")
ylabel("Angular rate (rad/s)")
legend("w_x","w_y","w_z")
grid on
saveas(gcf, "Images/Detumble_n.png")

% Plot Bdot controller output
% figure;
% plot(t, Bdotx)
% hold on
% plot(t, Bdoty)
% hold on 
% plot(t, Bdotz)
% title("Bdot controller output")
% xlabel("Time")
% ylabel("Magnetic moment ()")
% legend("m_x","m_y","m_z")
% grid on

import mlreportgen.report.*
import mlreportgen.dom.*

% Check if document exists, append if it does
docName = 'DETUMBLE_FIGURES.docx';
if isfile(docName)
    rpt = Document(docName, 'docx');
else
    rpt = Document(docName, 'docx');
end

% Add simulation case number (Update this as needed)
simCaseNumber = 01; % Update this number for each simulation case
par0 = Paragraph(sprintf('SIMULATION CASE NUMBER: %d', simCaseNumber));
par0.Style = {Bold(true), FontSize('14pt')};
append(rpt, par0);

% Add description paragraph (Update this as needed)
par1 = Paragraph('Detumbling w/o magnetorquer. Initial conditions are: w = [5,1,2]. Tuning: K = 0.5,0.5,0.2');
append(rpt, par1);

% Add first figure
fig1 = Image('Images/Detumble_s.png');
fig1.Height = "13cm";
fig1.Width = "16cm";
append(rpt, fig1);

% Add description for second figure
par2 = Paragraph('Same plot but with orbital periods');
append(rpt, par2);

% Add second figure
fig2 = Image('Images/Detumble_n.png');
fig2.Height = "13cm";
fig2.Width = "16cm";
append(rpt, fig2);

% Add description for Bdot plot
% par3 = Paragraph('Plot for bdot output controller');
% append(rpt, par3);

% Close and save the document
close(rpt);

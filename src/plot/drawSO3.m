clc;close all;clear all;

r1 = [2.213543189415407;-2.232440572204550;0.010336490115222];
r2 = [2.237099332166530;-2.203636935670440;-0.015064533201445];
r3 = [-2.209845834987348;2.228267503535384;-0.014739710318185];

R(:,:,1) = expSO3(r1);
R(:,:,2) = expSO3(r2);
R(:,:,3) = expSO3(r3);

AxisLength = 1;
ShowArrowHead = false;
ShowLegend = true;
p = zeros(3,3);
R = R(:, :, 1:end) * AxisLength;
    
fig = figure('NumberTitle', 'off', 'Name', '6DOF Animation');
set(gca, 'drawmode', 'fast');
lighting phong;
set(gcf, 'Renderer', 'zbuffer');
hold on;
xlim([-1.5*AxisLength, 1.5*AxisLength])
axis equal;
grid minor;
title(' Calibration result: {\bf R}_{\rm B}^{\rm C}');
Xlabel = 'x';
Ylabel = 'y';
Zlabel = 'z';
xlabel(Xlabel,'Interpreter','latex');
ylabel(Ylabel,'Interpreter','latex');
zlabel(Zlabel,'Interpreter','latex');

x(1) = p(1,1);
y(1) = p(1,2);
z(1) = p(1,3);
ox(1) = x(1);
oy(1) = y(1);
oz(1) = z(1);
ux(1) = R(1,1,1:1);
vx(1) = R(2,1,1:1);
wx(1) = R(3,1,1:1);
uy(1) = R(1,2,1:1);
vy(1) = R(2,2,1:1);
wy(1) = R(3,2,1:1);
uz(1) = R(1,3,1:1);
vz(1) = R(2,3,1:1);
wz(1) = R(3,3,1:1);
    
% Create graphics handles
orgHandle = plot3(x, y, z, 'k.');
if(ShowArrowHead)
    ShowArrowHeadStr = 'on';
else
    ShowArrowHeadStr = 'off';
end
    
% Create legend
if(ShowLegend)
%     legend('Origin', 'X', 'Y', 'Z');
end
    
Spin = 120;
View = [(100:(Spin/(3-1)):(100+Spin))', 10*ones(3, 1)]

cmap = lines(3);

for i = 1:size(R,3)
    odx = p(i,1);
    ody = p(i,2);
    odz = p(i,3);
    udx(1) = R(1,1,i);
    vdx(1) = R(2,1,i);
    wdx(1) = R(3,1,i);
    udy(1) = R(1,2,i);
    vdy(1) = R(2,2,i);
    wdy(1) = R(3,2,i);
    udz(1) = R(1,3,i);
    vdz(1) = R(2,3,i);
    wdz(1) = R(3,3,i);
    quivXhandle{i} = quiver3(odx, ody, odz, udx, vdx, wdx,  'Color', cmap(i,:), 'ShowArrowHead', 'on',  ...
                         'MaxHeadSize', 0.3, 'AutoScale', 'off','LineWidth',2);
    quivYhandle{i} = quiver3(odx, ody, odz, udy, vdy, wdy,  'Color', cmap(i,:), 'ShowArrowHead', 'on',  ...
                         'MaxHeadSize', 0.3, 'AutoScale', 'off','LineWidth',2);
    quivZhandle{i} = quiver3(odx, ody, odz, udz, vdz, wdz,  'Color', cmap(i,:), 'ShowArrowHead', 'on',  ...
                         'MaxHeadSize', 0.3, 'AutoScale', 'off','LineWidth',2);
end
% quivXhandle = quiver3(ox, oy, oz, ux, vx, wx,  'r', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off','LineWidth',1);
% quivYhandle = quiver3(ox, oy, oz, uy, vy, wy,  'g', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off','LineWidth',1);
% quivZhandle = quiver3(ox, ox, oz, uz, vz, wz,  'b', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off','LineWidth',1);
view(View(1, :));
if(ShowLegend)
    legend([quivXhandle{1},quivXhandle{2},quivXhandle{3}],{'Jia13','Jia14','Proposed'},'Interpreter','latex','Location','best');
end
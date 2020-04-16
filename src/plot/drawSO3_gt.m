clc;close all;clear all;



r3 = [-0.002058415098345;-0.001902305847459;-1.569983162699367];
r2 = [0.000008983662884;-0.000031304192247;-1.570487025713360];
r1 = [0.015281512148864;0.016400431748219;-1.571440605526156];

R(:,:,1) = expSO3(r1);
R(:,:,2) = [-0.0712 0.9974 -0.0130; -0.9954 -0.0719 -0.0640; -0.0648 0.0084 0.9979];

R(:,:,3) = expSO3(r2);
R(:,:,4) = [ 0.0288 0.9996 -0.0011; -0.9996 0.0288 -0.0030;-0.0029 0.0012 1.0000];

R(:,:,5) = expSO3(r3);
R(:,:,6) = [ 0.0213 0.9998 0.0024; -0.9996 0.0212 0.0164; 0.0164 -0.0028 0.9999];

%;%, ...
%                      [ 0.0288 0.9996 -0.0011; -0.9996 0.0288 -0.0030;-0.0029 0.0012 1.0000], ...
%                      [ 0.0213 0.9998 0.0024; -0.9996 0.0212 0.0164; 0.0164 -0.0028 0.9999]};

AxisLength = 1;
ShowArrowHead = false;
ShowLegend = true;
p = zeros(6,3);
R = R(:, :, 1:end) * AxisLength;
    
fig = figure('NumberTitle', 'off', 'Name', '6DOF Animation');
set(gca, 'drawmode', 'fast');
lighting phong;
set(gcf, 'Renderer', 'zbuffer');

title('$\mathbf{R}_{I}^{C}$','Interpreter','latex');
Xlabel = 'X';
Ylabel = 'Y';
Zlabel = 'Z';
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
% orgHandle = plot3(x, y, z, 'k.');
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
View = [(100:(Spin/(4-1)):(100+Spin))', 10*ones(4, 1)]

cmap = lines(4);

k = 1;
titlename = {'rccar','rotation','walk'};
for i = 1:2:size(R,3)
    subplot(1,3,k);
    
    xlim([-1.5*AxisLength, 1.5*AxisLength])
    axis equal;
    grid minor;
    
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
    quivXhandle{i} = quiver3(odx, ody, odz, udx, vdx, wdx,  'r', 'ShowArrowHead', 'on',  ...
                         'MaxHeadSize', 0.3, 'AutoScale', 'off','LineWidth',2);hold on;

    quivYhandle{i} = quiver3(odx, ody, odz, udy, vdy, wdy,  'r', 'ShowArrowHead', 'on',  ...
                         'MaxHeadSize', 0.3, 'AutoScale', 'off','LineWidth',2);
    quivZhandle{i} = quiver3(odx, ody, odz, udz, vdz, wdz,  'r', 'ShowArrowHead', 'on',  ...
                         'MaxHeadSize', 0.3, 'AutoScale', 'off','LineWidth',2);
    odx = p(i+1,1);
    ody = p(i+12);
    odz = p(i+1,3);
    udx(1) = R(1,1,i+1);
    vdx(1) = R(2,1,i+1);
    wdx(1) = R(3,1,i+1);
    udy(1) = R(1,2,i+1);
    vdy(1) = R(2,2,i+1);
    wdy(1) = R(3,2,i+1);
    udz(1) = R(1,3,i+1);
    vdz(1) = R(2,3,i+1);
    wdz(1) = R(3,3,i+1);
    quivXhandle{i+1} = quiver3(odx, ody, odz, udx, vdx, wdx,  'b', 'ShowArrowHead', 'on',  ...
                         'MaxHeadSize', 0.3, 'AutoScale', 'off','LineWidth',2);hold on;

    quivYhandle{i+1} = quiver3(odx, ody, odz, udy, vdy, wdy,  'b', 'ShowArrowHead', 'on',  ...
                         'MaxHeadSize', 0.3, 'AutoScale', 'off','LineWidth',2);
    quivZhandle{i+1} = quiver3(odx, ody, odz, udz, vdz, wdz,  'b', 'ShowArrowHead', 'on',  ...
                         'MaxHeadSize', 0.3, 'AutoScale', 'off','LineWidth',2);
                     
    if k == 1
        legend([quivXhandle{2},quivXhandle{1}],{'reference','estimation'},'Interpreter','latex');
    end
    title(titlename{k},'Interpreter','latex');
    k = k + 1;
    
    disp(norm(logSO3(R(:,:,i+1)'*R(:,:,i))))
end
% quivXhandle = quiver3(ox, oy, oz, ux, vx, wx,  'r', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off','LineWidth',1);
% quivYhandle = quiver3(ox, oy, oz, uy, vy, wy,  'g', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off','LineWidth',1);
% quivZhandle = quiver3(ox, ox, oz, uz, vz, wz,  'b', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off','LineWidth',1);
view(View(1, :));
% if(ShowLegend)
%     legend([quivXhandle{4},quivXhandle{1},quivXhandle{2},quivXhandle{3}],{'truth','rccar','rotation','walk'},'Interpreter','latex');
% end
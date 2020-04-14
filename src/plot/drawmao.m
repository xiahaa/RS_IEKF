clc;clear;close all;
syms x
y=ft(x);
ezplot(y);
grid minor
title('Mapping: 0.5+f\_{sig}(x)','Interpreter','latex');
ylabel('y','Interpreter','latex')
xlabel('x','Interpreter','latex')
hold on;
plot([0,0],[0,ft(0)],'r--','LineWidth',2);
xranges = xlim;
plot([xranges(1),0],[ft(0),ft(0)],'r--','LineWidth',2);
plot(0,ft(0),'rx','Markersize',20);


function dy = dft(x)
    dy = exp(x)/(exp(x) + 1) - exp(2*x)/(exp(x) + 1)^2;
end

function y = ft(x)
    fsig = @(x) (exp(x)/(exp(x)+1));
    y = (0.5+fsig(x));
end
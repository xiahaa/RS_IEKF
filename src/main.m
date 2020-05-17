%% demo code
clear all;clc;close all;

addpath(genpath('../thirdparty/gyro_camera_calibration/'));
addpath(genpath('./'));

%% call thirdparty program to load their data.
% this part is the original code provided by Jia Chao.
% Jia, Chao, and Brian L. Evans. 
% "Online calibration and synchronization of cellphone camera and gyroscope." 
% 2013 IEEE Global Conference on Signal and Information Processing. IEEE, 2013.
% 
[gyrostamp, gyrogap, anglev, firststamp] = readgyro('data/gyro.txt');
[framestamp, framegap] = readts('data/framestamps.txt');
framestart = 1;
frameend = length(framestamp);
% read matched features
[match_idx, match_x1, match_x2] = read_feature('data');

%% change initial value or filter type here.
para = config;

endidx = frameend-10;
warning off;

[mean_est, var_est, npara] = para.flt(match_idx, match_x1, match_x2, gyrostamp, gyrogap, anglev, framestamp, para, endidx);

% the calibrated parameters are saved in "npara"
npara

%% draw
t_span = framestamp(1:size(mean_est,1));
t_span = t_span(1:2:end);
figure(1);
k = 1;
subplot(2,2,k);
for i = 1:2
	sigma = sqrt(var_est(1:2:end,i));
	v_low  = mean_est(1:2:end,i) - sigma;
	v_high = mean_est(1:2:end,i) + sigma;
	fill([t_span;t_span(end:-1:1)],[v_low;v_high(end:-1:1)],[0 116 186]/255,'facealpha',.3);
	hold on
	plot(t_span,mean_est(1:2:end,i),'Color','r','LineWidth',1.5);hold on;grid on;
end
if k == 1
	%legend({'1 sigma envelope','posterior'},'Interpreter','latex');
	title('$(c_u, c_v)$','Interpreter','latex');
	xlabel('t: (s)','Interpreter','latex');
	ylabel('pixel','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex')
end
set(gca,'TickLabelInterpreter','latex')

name={'$t_r$', '$t_d$', '$f$'};
ylbs={'(s)','(s)','pixel'};
for i = 3:5
	k = k+1;
    subplot(2,2,k);
    sigma = sqrt(var_est(1:2:end,i));
	v_low  = mean_est(1:2:end,i) - sigma;
	v_high = mean_est(1:2:end,i) + sigma;
	fill([t_span;t_span(end:-1:1)],[v_low;v_high(end:-1:1)],[0 116 186]/255,'facealpha',.3);
	hold on
    plot(t_span,mean_est(1:2:end,i),'Color','r','LineWidth',1.5);hold on;grid on;
    title(name{i-2},'Interpreter','latex');
    xlabel('t: (s)','Interpreter','latex');
	ylabel(ylbs{i-2},'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    if k == 2
        legend({'1 sigma envelope','posterior'},'Interpreter','latex');
    end
end
warning on;

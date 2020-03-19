% Online calibration & synchronization of RS camera and gyro
% Calibration is done using an EKF with coplanarity constraints as measurements

clear all;clc;close all;

addpath(genpath('./'));
addpath(genpath('../thirdparty/gyro_camera_calibration/'));

para = config_dist(2);
[gyrostamp, gyrogap, anglev, firststamp] = readgyro_new(fullfile(para.basepath, 'gyro.txt'));
[framestamp, framegap] = readts_new(fullfile(para.basepath,'framestamps.txt'));

framestart = 1;

% read matched features
[match_idx, match_x1, match_x2] = read_feature(para.basepath);
warning off;
frameend = match_idx(end) - 2;

endidx = frameend-10;

% try
%     [mean_est1, var_est1, npara1] = ekf_epipolar(match_idx, match_x1, match_x2, gyrostamp, gyrogap, anglev, framestamp, para, endidx);
% catch
%     warning('ekf_epipolar fail!');
% end
% try
%     [mean_est2, var_est2, npara2] = ekf_epipolar_analytic(match_idx, match_x1, match_x2, gyrostamp, gyrogap, anglev, framestamp, para, endidx);
% catch
%     warning('ekf_epipolar_analytic fail!');
% end
if isfield(para,'fix')
    [mean_est, var_est, npara] = para.flt(match_idx, match_x1, match_x2, gyrostamp, gyrogap, anglev, framestamp, para, endidx,1);
else
    [mean_est, var_est, npara] = para.flt(match_idx, match_x1, match_x2, gyrostamp, gyrogap, anglev, framestamp, para, endidx);
end
% the calibrated parameters are saved in "npara"
npara

%% draw
t_span = framestamp(1:size(mean_est,1));
t_span = t_span(1:2:end);

if ~isfield(para,'fix')
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
    xlim([t_span(1) t_span(end)+0.5]);
    if k == 1
        legend({'1 sigma envelope','posterior'});
        title('(c_u, c_v)');
        xlabel('t: (s)');
        ylabel('pixel');
    end

    name={'t_r', 't_d', 'f'};
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
        title(name{i-2});
        xlabel('t: (s)');
        ylabel(ylbs{i-2});
        xlim([t_span(1) t_span(end)+0.5]);
    end
else
    figure(1);
    k = 1;
    name={'t_r', 't_d'};
    ylbs={'(s)','(s)'};
    for i = 1:2
        subplot(1,2,k);        k = k+1;
        sigma = sqrt(var_est(1:2:end,i));
        v_low  = mean_est(1:2:end,i) - sigma;
        v_high = mean_est(1:2:end,i) + sigma;
        fill([t_span;t_span(end:-1:1)],[v_low;v_high(end:-1:1)],[0 116 186]/255,'facealpha',.3);
        hold on
        plot(t_span,mean_est(1:2:end,i),'Color','r','LineWidth',1.5);hold on;grid on;
        title(name{i});
        xlabel('t: (s)');
        ylabel(ylbs{i});
        xlim([t_span(1) t_span(end)+0.5]);
    end
    subplot(1,2,1);
    plot([t_span(1), t_span(end)],[para.tr_th, para.tr_th],'k--','LineWidth',1.2);
    legend({'1 sigma envelope','posterior','truth'});
    subplot(1,2,2);
    plot([t_span(1), t_span(end)],[para.td_th, para.td_th],'k--','LineWidth',1.2);
    
    figure(2);
    Ric_th = para.ric_th;
    for i = 1:size(mean_est,1)
        Ric_est = expSO3(mean_est(i,6:8)');
        errs(i) = norm(logSO3(Ric_th'*Ric_est));
    end
    plot(framestamp(1:size(mean_est,1)),errs,'Color','r','LineWidth',1.5);hold on;grid minor;
    xlabel('t: (s)');
    ylabel('$\|\log(\mathbf{R}_{th}^T\bar{\mathbf{R}})\|$','Interpreter','latex');
end
warning on;

function [timestamp, timegap, anglev, firststamp] = readgyro_new( filepath )
    % read gyro information with timestamps
    [timestamp, anglevx, anglevy, anglevz] = textread(filepath, '%f,%f,%f,%f,');
    anglev = [anglevx, anglevy, anglevz];
    firststamp = timestamp(1);
    timestamp = timestamp - timestamp(1); % make the first timestamp as 0
    timegap = diff(timestamp);
    des_freq = 100;
    if mean(timegap) < (1/des_freq)
        % interpolation and descresing the frequency
        timestamp_new = timestamp(1):1/des_freq:timestamp(end);
        anglevx_new = interp1(timestamp,anglevx,timestamp_new,'spline');
        anglevy_new = interp1(timestamp,anglevy,timestamp_new,'spline');
        anglevz_new = interp1(timestamp,anglevz,timestamp_new,'spline');
        timestamp = timestamp_new;
        anglev = [anglevx_new', anglevy_new', anglevz_new'];
        timegap = diff(timestamp);
    end
end

function [timestamp, timegap] = readts_new(filepath)
    % read frame timestamps for each frame
    fid = fopen(filepath);
    timestamp = fscanf(fid,'%f');
    timegap = diff(timestamp);
    fclose(fid);
end

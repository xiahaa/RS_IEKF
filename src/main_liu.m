% Online calibration & synchronization of RS camera and gyro
% Calibration is done using an EKF with coplanarity constraints as measurements

clear all;clc;close all;

addpath(genpath('./'));
addpath(genpath('../thirdparty/gyro_camera_calibration/'));

para = config(2);
[gyrostamp, gyrogap, anglev, firststamp] = readgyro_new(fullfile(para.basepath, 'gyro.txt'));
[framestamp, framegap] = readts_new(fullfile(para.basepath,'framestamps.txt'));

framestart = 1;

% read matched features
[match_idx, match_x1, match_x2] = read_feature(para.basepath);
warning off;
frameend = match_idx(end) - 2;

endidx = frameend-10;
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

% Online calibration & synchronization of RS camera and gyro
% Calibration is done using an EKF with coplanarity constraints as measurements

clear all;
clc;
close all;

addpath(genpath('../thirdparty/gyro_camera_calibration/'));
addpath(genpath('../'));

[gyrostamp, gyrogap, anglev, firststamp] = readgyro('data/gyro.txt');
[framestamp, framegap] = readts('data/framestamps.txt');

framestart = 1;
frameend = length(framestamp);

% read matched features
[match_idx, match_x1, match_x2] = read_feature('data');

% set the initial values of parameters
para.ts = 0.025; % readout time of rolling shutter camera, 0.025 
para.td = 1.175;%1.175; % timestamp delay	
para.wd = [0 0 0]; % gyroscope bias
para.f = 700; % camera focal length
para.fx = 671; % camera focal length
para.fy = 671; % camera focal length
para.w = 720; % frame width (this is a fixed value)
para.h = 480; % frame height (this is a fixed value)
para.cx = (para.w-1)/2; % principal point
para.cy = (para.h-1)/2; % principal point
para.rcam = [pi/sqrt(2); -pi/sqrt(2); 0]; % relative orientation between camera and gyro
para.sigma = 1e-6; % the noise variance of gyro readings
para.pn = 1; % noise variance of feature detection

sigma_draw = [20,20,20,0.1,0.01];

N = 50;
failures1 = zeros(1,N);
failures2 = zeros(1,N);
failures3 = zeros(1,N);
res1 = zeros(N,5);
res2 = zeros(N,5);
res3 = zeros(N,5);

nomial.f = para.f;
nomial.cx = para.cx;
nomial.cy = para.cy;
nomial.ts = para.ts;
nomial.td = para.td;

for i = 1:N
    disp(i);
    draw = random_draw(sigma_draw);
    para.f = nomial.f + draw(1);
    para.cx = nomial.cx + draw(2);
    para.cy = nomial.cy + draw(3);
    para.ts = nomial.ts + draw(5);
    para.td = nomial.td + draw(4);

    % ekf filtering using transfer error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    endidx = frameend-10;
    try
        [mean_est1, var_est1, npara1] = ekf_epipolar(match_idx, match_x1, match_x2, gyrostamp, gyrogap, anglev, framestamp, para, endidx);
        res1(i,:) = [npara1.f,npara1.cx,npara1.cy,npara1.ts,npara1.td];
    catch
        failures1(i) = 1;
    end
    try
        [mean_est2, var_est2, npara2] = ekf_epipolar_analytic(match_idx, match_x1, match_x2, gyrostamp, gyrogap, anglev, framestamp, para, endidx);
        res2(i,:) = [npara2.f,npara2.cx,npara2.cy,npara2.ts,npara2.td];
    catch
        failures2(i) = 1;
    end
    % [mean_est, var_est, npara] = mhe_epipolar_analytic(match_idx, match_x1, match_x2, gyrostamp, gyrogap, anglev, framestamp, para, endidx);
    try
        [mean_est3, var_est3, npara3] = InEKF_rs_calib_rep(match_idx, match_x1, match_x2, gyrostamp, gyrogap, anglev, framestamp, para, endidx);
        res3(i,:) = [npara3.f,npara3.cx,npara3.cy,npara3.ts,npara3.td];
    catch
        failures3(i) = 1;
    end
    % [mean_est3, var_est3, npara3] = InEKF_rs_calib_fix_intrinsic(match_idx, match_x1, match_x2, gyrostamp, gyrogap, anglev, framestamp, para, endidx);
    % [mean_est3, var_est3, npara3] = InEKF_rs_calib_rep_exp(match_idx, match_x1, match_x2, gyrostamp, gyrogap, anglev, framestamp, para, endidx);
end
save('draw_res.mat','res1','res2','res3','failures1','failures2','failures3');

function rd = random_draw(sigma)
    rd = randn(1,length(sigma)).*sigma;
end
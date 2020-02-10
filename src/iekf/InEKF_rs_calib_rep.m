function [mean_est, var_est, npara] = InEKF_rs_calib_rep(match_idx, match_x1, match_x2, ...
    gyrostamp, gyrogap, anglev, framestamp, para, endidx)
%FUNCTION_NAME - runs InEKF for rolling shutter camera calibration.
%
% Syntax:  inde = update_inde(inde, x_k_k, p_k_k)
%
% Inputs:
%    match_idx - feature matching index.
%    match_x1 - feature matching position of frame i-1.
%    match_x2 - feature matching position of frame i.
%    gyrostamp - IMU timestamp.
%    gyrogap - gap time of two IMU samples.
%    anglev - angular velocity.
%    framestamp - img frame timestamp.
%    para - auxiliary parameter.
%    endidx - ending index.
%
% Outputs:
%    mean_est - logged mean values.
%    var_est - logged variance.
%    npara - final calibration parameters.
%
%
% Author: Xiao Hu
% Technical University of Denmark
% email: xiahaa@space.dtu.dk
% Jan 2020; Last revision: 31-Jan-2020
%------------- BEGIN CODE --------------
    h = para.h;
    w = para.w;
    gyro_sigma = para.sigma;% 
    pn = para.pn;% feature point noise
    std_ts = 0.01;
    std_td = 0.1;
    
    maxtd = std_td*3; % maximum range of parameter td is (+-) 0.03 s
    maxts = std_ts*3;
    % initial variance of different parameters
    oc_sigma = 6.67^2;
    ts_sigma = std_ts^2;
    td_sigma = std_td^2;
    f_sigma = 20^2;
    bias_sigma = 0.0067^2;
    rcam_sigma = 0.0067^2;
    
    Qw = gyro_sigma*eye(3); %% covariance matrix for angular velocity
    Qbw = 1e-6*eye(3);%% covariance matrix for bias of angular velocity
    Qs = blkdiag(Qw,Qbw);
    
    % initilize the EKF
    rcam = expSO3(para.rcam);
    x_k_k = [para.cx; para.cy; para.ts; para.td; para.f; para.wd'; rcam(:)];
    p_k_k = eye(11); 
    p_k_k(1,1) = oc_sigma; p_k_k(2,2) = oc_sigma; 
    p_k_k(3,3) = ts_sigma; p_k_k(4,4) = td_sigma; 
    p_k_k(5,5) = f_sigma; 
    p_k_k(6,6) = bias_sigma; p_k_k(7,7) = bias_sigma; p_k_k(8,8) = bias_sigma; 
    p_k_k(9,9) = rcam_sigma; p_k_k(10,10) = rcam_sigma; p_k_k(11,11) = rcam_sigma; 
    
    % group state
    inde.group_num = 1;
    inde.rot = 1:9; inde.cov_rot = 1:3;% latest index of rotation
    inde = update_inde(inde, x_k_k, p_k_k);
    
    R0 = eye(3);
    X(inde.group) = R0(:);
    X(inde.nongroup) = x_k_k';
    P = blkdiag(zeros(inde.cov_group(end),inde.cov_group(end)),p_k_k);
    
    frame_inde = 1;
    epipolar_num = [];
    ts_seq = [gyrostamp(1)];
    match_start = 1;
    upid = 1;
    
    %% filter start
    for ind = 1:length(gyrogap)
        tic;
        % time step
        dt = gyrogap(ind);
        % read gyro data and remove measured bias. 
        w=anglev(ind,:)-X(inde.wd);
        % fix handedness and orientation of gyro measuremet. 
        w = reshape(X(inde.rcam),3,3)*w';
        
        % prediction of current rotation matrix
        Rhat = reshape(X(inde.rot),3,3);
        Rpred = Rhat * expSO3(w*dt);
        % the rest state stay the same, so omit copy to save some time.
        
        % state augmentation
        X = [Rpred(:)' X];
        P = [P(inde.cov_rot,inde.cov_rot) P(inde.cov_rot,:);P(:,inde.cov_rot) P];
        inde.group_num = inde.group_num + 1;
        inde = update_inde(inde, x_k_k, p_k_k);
        ts_seq(end+1) = gyrostamp(ind+1);
        
        % covariance prediction
        % in this case, only the position term is nonzero
        F = eye(inde.cov_nongroup(end),inde.cov_nongroup(end));
        Ric = reshape(X(inde.rcam),3,3);
        F(inde.cov_rot(1):inde.cov_rot(end), inde.cov_wd(1):inde.cov_wd(end)) = -Rpred*Ric*dt;
        F(inde.cov_rot(1):inde.cov_rot(end), inde.cov_rcam(1):inde.cov_rcam(end)) = -Rpred*skewm(w)*dt;
        %         
        G = zeros(inde.cov_nongroup(end),size(Qs,1));
        Jrw = leftJ(w*dt);
        G(inde.cov_rot,1:3) = Rhat*Jrw*Ric;
        G(inde.cov_wd,4:6) = eye(3);
        % 
        P = F*P*F' + G*Qs*G'*dt;
        
        %% update part
        if isempty(epipolar_num)            
            lb_aligned_time = framestamp(frame_inde) + para.td-maxtd; 
            if ts_seq(end) < lb_aligned_time
                remov_can = find(ts_seq < lb_aligned_time);
                remov_can = remov_can(remov_can~=length(remov_can));% keep the latest
                ts_seq(remov_can) = [];
                remov_can = inde.group_num+1 - remov_can;
                canid = 1:inde.group_num;
                canid(remov_can) = [];
                %% delete useless states to keep the size of the filter
                X = [X(inde.group(1:canid(end)*9)) X(inde.nongroup)];
                P = [P(inde.cov_group(1:canid(end)*3),inde.cov_group(1:canid(end)*3)) P(inde.cov_group(1:canid(end)*3),inde.cov_nongroup(:));
                     P(inde.cov_nongroup(:), inde.cov_group(1:canid(end)*3)) P(inde.cov_nongroup(:), inde.cov_nongroup(:))];
                inde.group_num = inde.group_num - length(remov_can);
                inde.group_num = length(canid);
                inde = update_inde(inde, x_k_k, p_k_k);
            else
                epipolar_num(end+1) = frame_inde;
                frame_inde = frame_inde + 1;
            end
        else
            ub_aligned_time = framestamp(frame_inde) + para.td+maxtd+para.ts+maxts; 
            if ts_seq(end) > ub_aligned_time
                %% run ekf update
                epipolar_num(end+1) = frame_inde;
                frame_inde = frame_inde + 1;
                
                if epipolar_num(2) > endidx
                    break;
                end
                
                %% two frames ready, do epipolar constrained InEKF update
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % function ekf_info_z finds relevant features from matched result (a txt
                % file) and return two vectors, matchstart indicates the corresponding line
                % number in the txt file.
                [z_k, match_start, y_k] = ekf_info_z (epipolar_num(1), match_start, match_idx, match_x1, match_x2);
                % since each feature contains 2 position, so the feature number is the half
                % of the length of y_k
                fnum = length(y_k)/2;
                
                % EKF update
                if z_k ~= Inf
                    [mm, H_k, VRV_k] = rep_info_meas_analytical (X, inde, ts_seq, framestamp(epipolar_num(1)), framestamp(epipolar_num(2)),...
                        z_k, y_k, h, pn, para);
                    % traditional KF update
                    [X, P] = iekf_update_analytical (X, inde, P, H_k, VRV_k, zeros(length(mm),1), mm);
                else
                    % do nothing
                end
                
                % state deduction
                lb_aligned_time = framestamp(epipolar_num(2)) + para.td-maxtd;        
                remov_can = find(ts_seq < lb_aligned_time);
                remov_can = remov_can(remov_can~=length(ts_seq));% keep the latest
                ts_seq(remov_can) = [];
                remov_can = inde.group_num+1 - remov_can;
                canid = 1:inde.group_num;
                canid(remov_can) = [];
                %% delete useless states to keep the size of the filter
                X = [X(inde.group(1:canid(end)*9)) X(inde.nongroup)];
                P = [P(inde.cov_group(1:canid(end)*3),inde.cov_group(1:canid(end)*3)) P(inde.cov_group(1:canid(end)*3),inde.cov_nongroup(:));
                     P(inde.cov_nongroup(:), inde.cov_group(1:canid(end)*3)) P(inde.cov_nongroup(:), inde.cov_nongroup(:))];
                inde.group_num = length(canid);
                inde = update_inde(inde, x_k_k, p_k_k);
                
                % log
                var_est(upid,:) = [diag(P(inde.cov_nongroup,inde.cov_nongroup))'];
                x_k_display = [X(inde.nongroup(1:end-9)), logSO3(reshape(X(inde.nongroup(end-8:end)),3,3))'];
                mean_est(upid,:) = [X(inde.nongroup(1:end-9)), logSO3(reshape(X(inde.nongroup(end-8:end)),3,3))'];% convert rotation matrix to so3
                fprintf('idx = %d, ts = %f, td = %f, ocw = %f, och = %f, f = %f bias = %f %f %f rcamx = %f %f %f \n',epipolar_num(1), ...
                    x_k_display(3),x_k_display(4),x_k_display(1),x_k_display(2),x_k_display(5),x_k_display(6),x_k_display(7),x_k_display(8),...
                    x_k_display(9), x_k_display(10), x_k_display(11));        
                upid = upid + 1;
                % state deduction will be done in next frame.
                epipolar_num = [];
            else
                %% otherwise, wait
            end
        end
        runtime(ind) = toc;
    end
    x_k_display = [X(inde.nongroup(1:end-9)), logSO3(reshape(X(inde.nongroup(end-8:end)),3,3))'];
    npara.ts = x_k_display(3);
    npara.td = x_k_display(4);
    npara.wd = [x_k_display(6) x_k_display(7) x_k_display(8)];
    npara.f = x_k_display(5);
    npara.cx = x_k_display(1);
    npara.cy = x_k_display(2);
    npara.rcam = [x_k_display(9); x_k_display(10); x_k_display(11)];

    disp(['mean runtime is ', num2str(mean(runtime(runtime ~= 0)))]);
end






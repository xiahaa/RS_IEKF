function [mean_est, var_est, npara] = InEKF_rs_calib_rep_exp(match_idx, match_x1, match_x2, ...
	gyrostamp, gyrogap, anglev, framestamp, para, endidx, varargin)
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
    if isempty(varargin)
        fix_int = 0;
    else
        fix_int = 1;
    end
    
    use_pos_vel = 1;

    avg_freq = 1/mean(diff(framestamp));
    para.ts = para.ts;
    para.cx = para.cx;
    para.cy = para.cy;
    para.f = para.f;

    %% code section
    h = para.h;
    w = para.w;
    gyro_sigma = para.sigma;% 
    pn = para.pn;% feature point noise
    std_ts = 0.01;
    std_td = 0.1;
    
    maxtd = std_td*3; % maximum range of parameter td is (+-) 0.03 s
    maxts = std_ts*3;
    % initial variance of different parameters
    oc_sigma = 0.5^2;
    ts_sigma = std_ts^2;
    td_sigma = std_td^2;
    f_sigma = 0.5^2;
    bias_sigma = 0.0067^2;
    rcam_sigma = 0.0067^2;
    
    Qw = gyro_sigma*eye(3); %% covariance matrix for angular velocity
    Qbw = 1e-6*eye(3);%% covariance matrix for bias of angular velocity
    Qv = 1 * eye(3);%% covariance matrix for velocity noise    
    if use_pos_vel == 1
        Qs = blkdiag(Qw,Qbw,Qv);
    else
        Qs = blkdiag(Qw,Qbw);
    end
    % initilize the EKF
    rcam = expSO3(para.rcam);
    if fix_int == 0
        x_k_k = [0; 0; 0; 0; 0; para.wd'; rcam(:)];
        p_k_k = eye(11); 
        p_k_k(1,1) = oc_sigma; p_k_k(2,2) = oc_sigma; 
        p_k_k(3,3) = 0.2^2; p_k_k(4,4) = 0.1^2; 
        p_k_k(5,5) = f_sigma; 
        p_k_k(6,6) = bias_sigma; p_k_k(7,7) = bias_sigma; p_k_k(8,8) = bias_sigma; 
        p_k_k(9,9) = rcam_sigma; p_k_k(10,10) = rcam_sigma; p_k_k(11,11) = rcam_sigma; 
    else
        x_k_k = [0; 0; para.wd'; rcam(:)];
        p_k_k = eye(8); 
        p_k_k(1,1) = 0.4^2; p_k_k(2,2) = 0.4^2; 
        p_k_k(3,3) = bias_sigma; p_k_k(4,4) = bias_sigma; p_k_k(5,5) = bias_sigma; 
        p_k_k(6,6) = rcam_sigma; p_k_k(7,7) = rcam_sigma; p_k_k(8,8) = rcam_sigma; 
    end
        
    if isfield(para, 'dist')
        x_k_k = [x_k_k; para.dist'];
        dist_sigma = 0.03^2;
        p_k_k = blkdiag(p_k_k,dist_sigma.*eye(length(para.dist)));
    end
    
    % group state
    inde.group_num = 1;
    inde.rot = 1:9; inde.cov_rot = 1:3;% latest index of rotation
    if use_pos_vel == 1
        inde.pos = 10:12; inde.cov_pos = 4:6;% latest index of rotation
        inde.vel = 13:15; inde.cov_vel = 7:9;% latest index of rotation
    end
    if fix_int == 0
        inde.cxy = [];
    end
    inde = update_inde(inde, x_k_k, p_k_k,use_pos_vel);
    
    R0 = eye(3);
    if use_pos_vel == 1
        X(inde.group) = [R0(:);zeros(3,1);1*ones(3,1)];
        P = blkdiag(zeros(inde.cov_pos(end),inde.cov_pos(end)),1*eye(3),p_k_k);
    else
        X(inde.group) = R0(:);
        P = blkdiag(zeros(inde.cov_group(end),inde.cov_group(end)),p_k_k);
    end
    X(inde.nongroup) = x_k_k';
    
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
        if use_pos_vel == 1
            Phat = X(inde.pos);
            Vhat = X(inde.vel);
            Ppred = Phat + X(inde.vel) * dt;
            Vpred = Vhat;
        end
        % the rest state stay the same, so omit copy to save some time.
        
        % state augmentation
        if use_pos_vel == 1
            X = [[Rpred(:)' Ppred Vpred] X];
            P = [P(1:inde.cov_vel(end),1:inde.cov_vel(end)) P(1:inde.cov_vel(end),:);P(:,1:inde.cov_vel(end)) P];
        else
            X = [Rpred(:)' X];
            P = [P(inde.cov_rot,inde.cov_rot) P(inde.cov_rot,:);P(:,inde.cov_rot) P];
        end
        inde.group_num = inde.group_num + 1;
        inde = update_inde(inde, x_k_k, p_k_k,use_pos_vel);
        ts_seq(end+1) = gyrostamp(ind+1);
        
        % covariance prediction
        % in this case, only the position term is nonzero
        F = eye(inde.cov_nongroup(end),inde.cov_nongroup(end));
        Ric = reshape(X(inde.rcam),3,3);
        F(inde.cov_rot(1):inde.cov_rot(end), inde.cov_wd(1):inde.cov_wd(end)) = -Rhat*Ric*dt;
        F(inde.cov_rot(1):inde.cov_rot(end), inde.cov_rcam(1):inde.cov_rcam(end)) = -Rhat*skewm(w)*dt;
        if use_pos_vel == 1
            F(inde.cov_pos(1):inde.cov_pos(end),inde.cov_vel(1):inde.cov_vel(end)) = eye(3)*dt;

            F(inde.cov_vel(1):inde.cov_vel(end),inde.cov_wd(1):inde.cov_wd(end)) = -skewm(Vhat)*Rhat*Ric*dt;
            F(inde.cov_pos(1):inde.cov_pos(end),inde.cov_wd(1):inde.cov_wd(end)) = -skewm(Phat)*Rhat*Ric*dt;

            F(inde.cov_vel(1):inde.cov_vel(end),inde.cov_rcam(1):inde.cov_rcam(end)) = -skewm(Vhat)*Rhat*skewm(w)*dt;
            F(inde.cov_pos(1):inde.cov_pos(end),inde.cov_rcam(1):inde.cov_rcam(end)) = -skewm(Phat)*Rhat*skewm(w)*dt;
        end
        
        %         
        G = zeros(inde.cov_nongroup(end),size(Qs,1));
        Jrw = leftJ(w*dt);
        G(inde.cov_rot,1:3) = Rhat*Jrw*Ric;
        G(inde.cov_wd,4:6) = eye(3);
        G(inde.cov_rot,4:6) = Rhat*Jrw*Ric*dt;
        if use_pos_vel == 1
            G(inde.cov_vel,[1:3,7:9]) = [skewm(Vhat)*Rhat*Jrw*Ric, eye(3)];
            G(inde.cov_pos,[1:3,7:9]) = [skewm(Ppred)*Rhat*Jrw*Ric, eye(3)*dt];
            G(inde.cov_vel,4:6) = skewm(Vhat)*Rhat*Jrw*Ric*dt;
            G(inde.cov_pos,4:6) = skewm(Ppred)*Rhat*Jrw*Ric*dt;
        end
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
                if use_pos_vel == 1
                    X = [X(inde.group(1:canid(end)*15)) X(inde.nongroup)];
                    P = [P(inde.cov_group(1:canid(end)*9),inde.cov_group(1:canid(end)*9)) P(inde.cov_group(1:canid(end)*9),inde.cov_nongroup(:));
                         P(inde.cov_nongroup(:), inde.cov_group(1:canid(end)*9)) P(inde.cov_nongroup(:), inde.cov_nongroup(:))];
                else
                    X = [X(inde.group(1:canid(end)*9)) X(inde.nongroup)];
                    P = [P(inde.cov_group(1:canid(end)*3),inde.cov_group(1:canid(end)*3)) P(inde.cov_group(1:canid(end)*3),inde.cov_nongroup(:));
                         P(inde.cov_nongroup(:), inde.cov_group(1:canid(end)*3)) P(inde.cov_nongroup(:), inde.cov_nongroup(:))];
                end
                
                inde.group_num = inde.group_num - length(remov_can);
                inde.group_num = length(canid);
                inde = update_inde(inde, x_k_k, p_k_k,use_pos_vel);
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
                    if use_pos_vel == 1
                        [mm, H_k, VRV_k] = rep_info_meas_analytical2 (X, inde, ts_seq, framestamp(epipolar_num(1)), framestamp(epipolar_num(2)),...
                        z_k, y_k, h, pn, para);
                    else
                        [mm, H_k, VRV_k] = rep_info_meas_analytical (X, inde, ts_seq, framestamp(epipolar_num(1)), framestamp(epipolar_num(2)),...
                            z_k, y_k, h, pn, para);
                    end
                    % traditional KF update
                    [X, P] = iekf_update_analytical (X, inde, P, H_k, VRV_k, zeros(length(mm),1), mm, use_pos_vel);
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
                if use_pos_vel == 1
                    X = [X(inde.group(1:canid(end)*15)) X(inde.nongroup)];
                    P = [P(inde.cov_group(1:canid(end)*9),inde.cov_group(1:canid(end)*9)) P(inde.cov_group(1:canid(end)*9),inde.cov_nongroup(:));
                         P(inde.cov_nongroup(:), inde.cov_group(1:canid(end)*9)) P(inde.cov_nongroup(:), inde.cov_nongroup(:))];
                else
                    X = [X(inde.group(1:canid(end)*9)) X(inde.nongroup)];
                    P = [P(inde.cov_group(1:canid(end)*3),inde.cov_group(1:canid(end)*3)) P(inde.cov_group(1:canid(end)*3),inde.cov_nongroup(:));
                         P(inde.cov_nongroup(:), inde.cov_group(1:canid(end)*3)) P(inde.cov_nongroup(:), inde.cov_nongroup(:))];
                end
                inde.group_num = length(canid);
                inde = update_inde(inde, x_k_k, p_k_k,use_pos_vel);
                
                % log
                if isfield(inde, 'cxy')
                    ccx = para.cx*dft(X(inde.cxy(1)));
                    ccy = para.cy*dft(X(inde.cxy(2)));
                    ccf = para.f*dft(X(inde.f));
                end
                ccc = para.ts*dft(X(inde.ts));
                ccd = para.td*dft(X(inde.td));
                var_est(upid,:) = [diag(P(inde.cov_nongroup,inde.cov_nongroup))'];
                if isfield(inde, 'cxy')
                    var_est(upid,1) = ccx * var_est(upid,1) * ccx';
                    var_est(upid,2) = ccy * var_est(upid,2) * ccy';
                    var_est(upid,5) = ccf * var_est(upid,5) * ccf';
                    var_est(upid,3) = ccc * var_est(upid,3) * ccc';
                    var_est(upid,4) = ccd * var_est(upid,4) * ccd';
                else
                    var_est(upid,1) = ccc * var_est(upid,1) * ccc';
                    var_est(upid,2) = ccd * var_est(upid,2) * ccd';
                end
                x_k_display = [X([inde.nongroup(1):inde.rcam(1)-1]), logSO3(reshape(X(inde.rcam),3,3))', X([inde.rcam(end)+1:inde.nongroup(end)])];
                if isfield(inde, 'cxy')
                    x_k_display(1) = para.cx * ft(X(inde.cxy(1)));
                    x_k_display(2) = para.cy * ft(X(inde.cxy(2)));
                    x_k_display(5) = para.f * ft(X(inde.f));
                    x_k_display(3) = para.ts * ft(X(inde.ts));
                    x_k_display(4) = para.td * ft(X(inde.td));
                    fprintf('idx = %d, ts = %f, td = %f, ocw = %f, och = %f, f = %f bias = %f %f %f rcamx = %f %f %f \n',epipolar_num(1), ...
                        x_k_display(3),x_k_display(4),x_k_display(1),x_k_display(2),x_k_display(5),x_k_display(6),x_k_display(7),x_k_display(8),...
                        x_k_display(9), x_k_display(10), x_k_display(11));        
                else
                    x_k_display(1) = para.ts * ft(X(inde.ts));
                    x_k_display(2) = para.td * ft(X(inde.td));
                    fprintf('idx = %d, ts = %f, td = %f, bias = %f %f %f rcamx = %f %f %f \n',epipolar_num(1), ...
                        x_k_display(1),x_k_display(2),x_k_display(3),x_k_display(4),x_k_display(5),...
                        x_k_display(6), x_k_display(7), x_k_display(8));        
                end
                mean_est(upid,:) = x_k_display;

                upid = upid + 1;
                % state deduction will be done in next frame.
                epipolar_num = [];
            else
                %% otherwise, wait
            end
        end
        runtime(ind) = toc;
    end
    x_k_display = [X([inde.nongroup(1):inde.rcam(1)-1]), logSO3(reshape(X(inde.rcam),3,3))', X([inde.rcam(end)+1:inde.nongroup(end)])];
    if isfield(inde, 'cxy')
        x_k_display(1) = para.cx * ft(X(inde.cxy(1)));
        x_k_display(2) = para.cy * ft(X(inde.cxy(2)));
        x_k_display(5) = para.f * ft(X(inde.f));
        x_k_display(3) = para.ts * ft(X(inde.ts));
        x_k_display(4) = para.td * ft(X(inde.td));
        npara.ts = x_k_display(3);
        npara.td = x_k_display(4);
        npara.wd = [x_k_display(6) x_k_display(7) x_k_display(8)];
        npara.f = x_k_display(5);
        npara.cx = x_k_display(1);
        npara.cy = x_k_display(2);
        npara.rcam = [x_k_display(9); x_k_display(10); x_k_display(11)];
        if length(x_k_display)>11
            npara.k1 = x_k_display(12);
        end
        if length(x_k_display)>12
            npara.k2 = x_k_display(13);
        end
        if length(x_k_display)>13
            npara.k3 = x_k_display(14);
        end
    else
        x_k_display(1) = para.ts * ft(X(inde.ts));
        x_k_display(2) = para.td * ft(X(inde.td));
        npara.ts = x_k_display(1);
        npara.td = x_k_display(2);
        npara.wd = [x_k_display(3) x_k_display(4) x_k_display(5)];
        npara.rcam = [x_k_display(6); x_k_display(7); x_k_display(8)];
        if length(x_k_display)>8
            npara.k1 = x_k_display(9);
        end
        if length(x_k_display)>9
            npara.k2 = x_k_display(10);
        end
        if length(x_k_display)>10
            npara.k3 = x_k_display(11);
        end
    end
    disp(['mean runtime is ', num2str(mean(runtime(runtime ~= 0)))]);
end

function dy = dft(x)
    dy = exp(x)/(exp(x) + 1) - exp(2*x)/(exp(x) + 1)^2;
end

function y = ft(x)
    fsig = @(x) (exp(x)./(exp(x)+1));
    y = (0.5+fsig(x));
end

function [yhat, H, JRJt] = rep_info_meas_analytical(X, inde, localstamp, fta, ftb, fz, fy, h, puv,para)
    x_nongroup = X(inde.nongroup);
    % split current states
	if isfield(inde,'cxy')
	    ocx = para.cx * ft(x_nongroup(1));
	    ocy = para.cy * ft(x_nongroup(2));
        tr = para.ts * ft(x_nongroup(3));
	    td = para.td * ft(x_nongroup(4));
	    f = para.f * ft(x_nongroup(5));
    else
        tr = para.ts * ft(x_nongroup(1));
    	td = para.td * ft(x_nongroup(2));
    end
    
    if isfield(inde,'k1')
        k1 = X(inde.k1);
    else
        k1 = 0;
    end
    
    if isfield(inde,'k2')
        k2 = X(inde.k2);
    else
        k2 = 0;
    end
    
    if isfield(inde,'k3')
        k3 = X(inde.k3);
    else
        k3 = 0;
    end

    hnum = length(fy)/2;
    H = zeros(hnum, inde.cov_nongroup(end));
    yhat = zeros(hnum,1);
    
    fta = fta + td;      % aligned frame time
    ftb = ftb + td; 
    ta = fta + fy(2:2:end)./h.*tr;% aligned feature time
    tb = ftb + fz(2:2:end)./h.*tr;
    
    % since for rs camera, each feature will have its own rotation matrix
    % because they may lay on different row. z means features in new frame, y
    % means features in old frame. 
    %% to d
	if isfield(inde, 'cxy')
		K = [f 0 ocx;0 f ocy;0 0 1];Kinv = inv(K);
    	fy1 = [reshape(fy, 2, length(fy)/2); f*ones(1, length(fy)/2)]; fy1 = fy1 - [ocx;ocy;0];
    	fz1 = [reshape(fz, 2, length(fz)/2); f*ones(1, length(fz)/2)]; fz1 = fz1 - [ocx;ocy;0];
        
        r1 = ((fy1(1,:) ./ fy1(3,:)).^2 + (fy1(2,:) ./ fy1(3,:)).^2);
        distcorr1 = 1 + k1*r1 + k2 * r1.^2 + k3 * r1.^3;
        dfy1 = [distcorr1 .* fy1(1:2,:);fy1(3,:)];
        
        r2 = ((fz1(1,:) ./ fz1(3,:)).^2 + (fz1(2,:) ./ fz1(3,:)).^2);
        distcorr2 = 1 + k1*r2 + k2 * r2.^2 + k3 * r2.^3;
        dfz1 = [distcorr2 .* fz1(1:2,:);fz1(3,:)];
    else
		K = [para.fx 0 para.cx;0 para.fy para.cy;0 0 1];Kinv = inv(K);
    	fy1 = [reshape(fy, 2, length(fy)/2); ones(1, length(fy)/2)]; fy1 = Kinv*fy1;
    	fz1 = [reshape(fz, 2, length(fz)/2); ones(1, length(fz)/2)]; fz1 = Kinv*fz1;
        r1 = ((fy1(1,:) ./ fy1(3,:)).^2 + (fy1(2,:) ./ fy1(3,:)).^2);
        distcorr1 = 1 + k1*r1 + k2 * r1.^2 + k3 * r1.^3;
        dfy1 = [distcorr1 .* fy1(1:2,:);fy1(3,:)];
        
        r2 = ((fz1(1,:) ./ fz1(3,:)).^2 + (fz1(2,:) ./ fz1(3,:)).^2);
        distcorr2 = 1 + k1*r2 + k2 * r2.^2 + k3 * r2.^3;
        dfz1 = [distcorr2 .* fz1(1:2,:);fz1(3,:)];
	end
    %% normalization to e
    nfy1 = dfy1 ./ vecnorm(dfy1);
    nfz1 = dfz1 ./ vecnorm(dfz1);
    
    param = struct;
    param.fy = reshape(fy, 2, length(fy)/2); param.fz = reshape(fz, 2, length(fz)/2);
    param.fy1 = fy1; param.fz1 = fz1;
    param.dfy1 = dfy1; param.dfz1 = dfz1;
    param.nfy1 = nfy1; param.nfz1 = nfz1;
    param.r1 = r1; param.r2 = r2;
    param.distcorr1=distcorr1; param.distcorr2=distcorr2;
    param.k1 = k1;param.k2 = k2;param.k3 = k3;
    
	Hf = zeros(hnum, 4);
    JRJt = zeros(hnum, hnum);
    delta = 0.01;
    usehuber = 1;
    for i = 1:hnum
        %% compute intermediate variables 
        ind_1 = i;
        
        %% compute output & jacobian
        [a1, b1, jac_a1_x, jac_a1_uv, jac_b1_x, jac_b1_uv] = cons_a(X, inde, localstamp, ta, tb, ind_1, h, param, Kinv, para);

        % output
        mm = 1-dot(a1,b1);
        if usehuber == 1
            [yhat(i), costd] = huber_cost(mm, delta);
        else
            yhat(i) = mm; costd = 1;
        end
        % jacobian
        cons1 = costd*(-b1');
        cons2 = costd*(-a1');
                
        H(i,:) = cons1*jac_a1_x + cons2 * jac_b1_x;
        
        Hf(i,:) = cons1*jac_a1_uv + cons2*jac_b1_uv;
        
        JRJt(i,i) = (cons1*jac_a1_uv)*(cons1*jac_a1_uv)' + (cons2*jac_b1_uv)*(cons2*jac_b1_uv)';
        JRJt(i,i) = JRJt(i,i).*puv;
    end
    %% left null basis
%     [U,~,~] = svd(Hf);
%     AA = U(:,size(Hf,2)+1:end);
%     H = AA'*H;
%     JRJt = puv*eye(size(H,1));
%     yhat = AA'*yhat;
end

function [a, b, jac_a_x, jac_a_uv, jac_b_x, jac_b_uv] = cons_a(X, inde, localstamp, ta, tb, ind_1, h, param, Kinv, para)
    fy = param.fy; fz = param.fz;
    fy1 = param.fy1; fz1 = param.fz1;
    dfy1 = param.dfy1; dfz1 = param.dfz1;
    nfy1 = param.nfy1; nfz1 = param.nfz1;
    r1 = param.r1; r2 = param.r2;
    distcorr1 = param.distcorr1; distcorr2 = param.distcorr2;
    k1 = param.k1;k2 = param.k2;k3 = param.k3;

    % call to get relative rotation
    [Ry, dty, iy, jy, R0y, R1y, xity, sy] = geodesic_interpolate_rot(X, inde, localstamp, ta(ind_1));
    [Rz, dtz, iz, jz, R0z, R1z, xitz, sz] = geodesic_interpolate_rot(X, inde, localstamp, tb(ind_1));
    
    a = Ry*nfy1(:,ind_1);
    b = Rz*nfz1(:,ind_1);
        
    xi1 = logSO3(R0y'*R1y);
    Jl1inv = leftJinv(xi1); Jr1inv = rightJinv(xi1);
    Jr2 = rightJ(xity);
    
    jac_R0 = Ry*skewm(nfy1(:,ind_1))*Jr2*Jl1inv*R0y'*dty-skewm(a);
    jac_R1 = -Ry*skewm(nfy1(:,ind_1))*Jr2*Jr1inv*R1y'*dty;
    
    tmp1 = -Ry*skewm(nfy1(:,ind_1))*xi1;
    
    %% jacobian
    jac_a_x = zeros(3,inde.cov_nongroup(end));
	if isfield(inde, 'cxy')
        ccx = para.cx*dft(X(inde.cxy(1)));
        ccy = para.cy*dft(X(inde.cxy(2)));
        ccf = para.f*dft(X(inde.f));
        
    	tmp2 = Ry*jac_f5(dfy1(:,ind_1));
        
        jac_cx = [-(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)^2/(fy1(3,ind_1)^2)-distcorr1(ind_1); ...
                  -(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*fy1(2,ind_1)/(fy1(3,ind_1)^2); ...
                  0];
        
        jac_cy = [-(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*fy1(2,ind_1)/(fy1(3,ind_1)^2); ...
                  -(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(2,ind_1)^2/(fy1(3,ind_1)^2)-distcorr1(ind_1); ...
                  0];
              
        jac_f = [-(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*(fy1(1,ind_1)/fy1(3,ind_1)^3+fy1(2,ind_1)/fy1(3,ind_1)^3); ...
                  -(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(2,ind_1)*(fy1(1,ind_1)/fy1(3,ind_1)^3+fy1(2,ind_1)/fy1(3,ind_1)^3); ...
                  1] + distcorr1(ind_1)*[0;0;1];
        
    	jac_a_x(:,[inde.cov_cxy, inde.cov_f]) = [tmp2*ccx*jac_cx tmp2*ccy*jac_cy tmp2*ccf*jac_f];
    
        if isfield(inde, 'k1')
            jac_a_x(:,[inde.cov_k1]) = tmp2*[r1(ind_1)*fy1(1,ind_1);r1(ind_1)*fy1(2,ind_1);0];
        end
        
        if isfield(inde, 'k2')
            jac_a_x(:,[inde.cov_k2]) = tmp2*[r1(ind_1)^2*fy1(1,ind_1);r1(ind_1)^2*fy1(2,ind_1);0];
        end
        
        if isfield(inde, 'k3')
            jac_a_x(:,[inde.cov_k3]) = tmp2*[r1(ind_1)^3*fy1(1,ind_1);r1(ind_1)^3*fy1(2,ind_1);0];
        end
    
        jac_u = [(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)^2/(fy1(3,ind_1)^2)+distcorr1(ind_1); ...
                  (2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*fy1(2,ind_1)/(fy1(3,ind_1)^2); ...
                  0];

        jac_v = [(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*fy1(2,ind_1)/(fy1(3,ind_1)^2); ...
                 (2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(2,ind_1)^2/(fy1(3,ind_1)^2)+distcorr1(ind_1); ...
                  0];
	else
    	tmp2 = Ry*jac_f5(dfy1(:,ind_1));
        if isfield(inde, 'k1')
            jac_a_x(:,[inde.cov_k1]) = tmp2*[r1(ind_1)*fy1(1,ind_1);r1(ind_1)*fy1(2,ind_1);0];
        end
        
        if isfield(inde, 'k2')
            jac_a_x(:,[inde.cov_k2]) = tmp2*[r1(ind_1)^2*fy1(1,ind_1);r1(ind_1)^2*fy1(2,ind_1);0];
        end
        
        if isfield(inde, 'k3')
            jac_a_x(:,[inde.cov_k3]) = tmp2*[r1(ind_1)^3*fy1(1,ind_1);r1(ind_1)^3*fy1(2,ind_1);0];
        end
        f1 = para.fx;
        f2 = para.fy;
        jac_u = [(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)^2/(f1)+distcorr1(ind_1)/f1; ...
                  (2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*fy1(2,ind_1)/(f1); ...
                  0];

        jac_v = [(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*fy1(2,ind_1)/(f2); ...
                 (2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(2,ind_1)^2/(f2)+distcorr1(ind_1)/f2; ...
                  0];
    end
    
    jac_a_uv(:,1:4) = [tmp2*jac_u tmp2*jac_v zeros(3,1) zeros(3,1)]+[zeros(3,1) tmp1*1/h*sy zeros(3,1) zeros(3,1)];

    v = fy(2,ind_1);
    ccc = para.ts*dft(X(inde.ts))*v/h;    
    ccd = para.td*dft(X(inde.td));    
    jac_a_x(:,[inde.cov_ts, inde.cov_td]) = [tmp1*sy*ccc tmp1*sy*ccd];
    
    jac_a_x(:,inde.cov_group(iy*3-2:iy*3)) = jac_R0;
    jac_a_x(:,inde.cov_group(jy*3-2:jy*3)) = jac_R1;
    
    %% b
    xi1 = logSO3(R0z'*R1z);
    Jl1inv = leftJinv(xi1); Jr1inv = rightJinv(xi1);
    Jr2 = rightJ(xitz);
    
    jac_R0 = Rz*skewm(nfz1(:,ind_1))*Jr2*Jl1inv*R0z'*dtz-skewm(b);
    jac_R1 = -Rz*skewm(nfz1(:,ind_1))*Jr2*Jr1inv*R1z'*dtz;
    
    tmp1 = -Rz*skewm(nfz1(:,ind_1))*xi1;
    
    jac_b_x = zeros(3,inde.cov_nongroup(end));
	if isfield(inde, 'cxy')
        ccx = para.cx*dft(X(inde.cxy(1)));
        ccy = para.cy*dft(X(inde.cxy(2)));
        ccf = para.f*dft(X(inde.f));
	    tmp2 = Rz*jac_f5(dfz1(:,ind_1));
        
        jac_cx = [-(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)^2/(fz1(3,ind_1)^2)-distcorr2(ind_1); ...
                  -(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*fz1(2,ind_1)/(fz1(3,ind_1)^2); ...
                  0];
        
        jac_cy = [-(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*fz1(2,ind_1)/(fz1(3,ind_1)^2); ...
                  -(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(2,ind_1)^2/(fz1(3,ind_1)^2)-distcorr2(ind_1); ...
                  0];
              
        jac_f = [-(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*(fz1(1,ind_1)/fz1(3,ind_1)^3+fz1(2,ind_1)/fz1(3,ind_1)^3); ...
                  -(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(2,ind_1)*(fz1(1,ind_1)/fz1(3,ind_1)^3+fz1(2,ind_1)/fz1(3,ind_1)^3); ...
                  1] + distcorr2(ind_1)*[0;0;1];
  
    	jac_b_x(:,[inde.cov_cxy, inde.cov_f]) = [tmp2*ccx*jac_cx tmp2*ccy*jac_cy tmp2*ccf*jac_f];
                                                   
        if isfield(inde, 'k1')
            jac_b_x(:,[inde.cov_k1]) = tmp2*[r2(ind_1)*fz1(1,ind_1);r2(ind_1)*fz1(2,ind_1);0];
        end
        
        if isfield(inde, 'k2')
            jac_b_x(:,[inde.cov_k2]) = tmp2*[r2(ind_1)^2*fz1(1,ind_1);r2(ind_1)^2*fz1(2,ind_1);0];
        end
        
        if isfield(inde, 'k3')
            jac_b_x(:,[inde.cov_k3]) = tmp2*[r2(ind_1)^3*fz1(1,ind_1);r2(ind_1)^3*fz1(2,ind_1);0];
        end
        
        jac_u = [(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)^2/(fz1(3,ind_1)^2)+distcorr2(ind_1); ...
                  (2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*fz1(2,ind_1)/(fz1(3,ind_1)^2); ...
                  0];

        jac_v = [(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*fz1(2,ind_1)/(fz1(3,ind_1)^2); ...
                 (2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(2,ind_1)^2/(fz1(3,ind_1)^2)+distcorr2(ind_1); ...
                  0];
    else
	    tmp2 = Rz*jac_f5(dfz1(:,ind_1));
        if isfield(inde, 'k1')
            jac_b_x(:,[inde.cov_k1]) = tmp2*[r2(ind_1);r2(ind_1);0];
        end
        
        if isfield(inde, 'k2')
            jac_b_x(:,[inde.cov_k2]) = tmp2*[r2(ind_1)^2;r2(ind_1)^2;0];
        end
        
        if isfield(inde, 'k3')
            jac_b_x(:,[inde.cov_k3]) = tmp2*[r2(ind_1)^3;r2(ind_1)^3;0];
        end
        
        f1 = para.fx;
        f2 = para.fy;
        jac_u = [(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)^2/(f1)+distcorr2(ind_1)/f1; ...
                  (2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*fz1(2,ind_1)/(f1); ...
                  0];

        jac_v = [(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*fz1(2,ind_1)/(f2); ...
                 (2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(2,ind_1)^2/(f2)+distcorr2(ind_1)/f2; ...
                  0];
    end
    
    jac_b_uv(:,1:4) = [zeros(3,1) zeros(3,1) tmp2*jac_u tmp2*jac_v] + [zeros(3,1) zeros(3,1) zeros(3,1) tmp1*1/h*sz];
    
    v = fz(2,ind_1);
    ccc = para.ts*dft(X(inde.ts))*v/h;   
    ccd = para.td*dft(X(inde.td));  
    jac_b_x(:,[inde.cov_ts, inde.cov_td]) = [tmp1*sz*ccc tmp1*sz*ccd];
    
    jac_b_x(:,inde.cov_group(iz*3-2:iz*3)) = jac_R0;
    jac_b_x(:,inde.cov_group(jz*3-2:jz*3)) = jac_R1;
end

function [yhat, H, JRJt] = rep_info_meas_analytical2(X, inde, localstamp, fta, ftb, fz, fy, h, puv,para)
    x_nongroup = X(inde.nongroup);
    % split current states
	if isfield(inde,'cxy')
	    ocx = para.cx * ft(x_nongroup(1));
	    ocy = para.cy * ft(x_nongroup(2));
        tr = para.ts * ft(x_nongroup(3));
	    td = para.td * ft(x_nongroup(4));
	    f = para.f * ft(x_nongroup(5));
    else
        tr = para.ts * ft(x_nongroup(1));
    	td = para.td * ft(x_nongroup(2));
    end
    
    if isfield(inde,'k1')
        k1 = X(inde.k1);
    else
        k1 = 0;
    end
    
    if isfield(inde,'k2')
        k2 = X(inde.k2);
    else
        k2 = 0;
    end
    
    if isfield(inde,'k3')
        k3 = X(inde.k3);
    else
        k3 = 0;
    end

    hnum = length(fy)/2;
    H = zeros(hnum, inde.cov_nongroup(end));
    yhat = zeros(hnum,1);
    
    fta = fta + td;      % aligned frame time
    ftb = ftb + td; 
    ta = fta + fy(2:2:end)./h.*tr;% aligned feature time
    tb = ftb + fz(2:2:end)./h.*tr;
    
    % since for rs camera, each feature will have its own rotation matrix
    % because they may lay on different row. z means features in new frame, y
    % means features in old frame. 
    %% to d
	if isfield(inde, 'cxy')
		K = [f 0 ocx;0 f ocy;0 0 1];Kinv = inv(K);
    	fy1 = [reshape(fy, 2, length(fy)/2); f*ones(1, length(fy)/2)]; fy1 = fy1 - [ocx;ocy;0];
    	fz1 = [reshape(fz, 2, length(fz)/2); f*ones(1, length(fz)/2)]; fz1 = fz1 - [ocx;ocy;0];
        
        r1 = ((fy1(1,:) ./ fy1(3,:)).^2 + (fy1(2,:) ./ fy1(3,:)).^2);
        distcorr1 = 1 + k1*r1 + k2 * r1.^2 + k3 * r1.^3;
        dfy1 = [distcorr1 .* fy1(1:2,:);fy1(3,:)];
        
        r2 = ((fz1(1,:) ./ fz1(3,:)).^2 + (fz1(2,:) ./ fz1(3,:)).^2);
        distcorr2 = 1 + k1*r2 + k2 * r2.^2 + k3 * r2.^3;
        dfz1 = [distcorr2 .* fz1(1:2,:);fz1(3,:)];
    else
		K = [para.fx 0 para.cx;0 para.fy para.cy;0 0 1];Kinv = inv(K);
    	fy1 = [reshape(fy, 2, length(fy)/2); ones(1, length(fy)/2)]; fy1 = Kinv*fy1;
    	fz1 = [reshape(fz, 2, length(fz)/2); ones(1, length(fz)/2)]; fz1 = Kinv*fz1;
        r1 = ((fy1(1,:) ./ fy1(3,:)).^2 + (fy1(2,:) ./ fy1(3,:)).^2);
        distcorr1 = 1 + k1*r1 + k2 * r1.^2 + k3 * r1.^3;
        dfy1 = [distcorr1 .* fy1(1:2,:);fy1(3,:)];
        
        r2 = ((fz1(1,:) ./ fz1(3,:)).^2 + (fz1(2,:) ./ fz1(3,:)).^2);
        distcorr2 = 1 + k1*r2 + k2 * r2.^2 + k3 * r2.^3;
        dfz1 = [distcorr2 .* fz1(1:2,:);fz1(3,:)];
	end
    %% normalization to e
    nfy1 = dfy1 ./ vecnorm(dfy1);
    nfz1 = dfz1 ./ vecnorm(dfz1);
    
    param = struct;
    param.fy = reshape(fy, 2, length(fy)/2); param.fz = reshape(fz, 2, length(fz)/2);
    param.fy1 = fy1; param.fz1 = fz1;
    param.dfy1 = dfy1; param.dfz1 = dfz1;
    param.nfy1 = nfy1; param.nfz1 = nfz1;
    param.r1 = r1; param.r2 = r2;
    param.distcorr1=distcorr1; param.distcorr2=distcorr2;
    param.k1 = k1;param.k2 = k2;param.k3 = k3;
    
	Hf = zeros(hnum, 4);
    JRJt = zeros(hnum, hnum);
    delta = 1;
    usehuber = 1;
    for i = 1:hnum
        %% compute intermediate variables 
        ind_1 = i;
        
        %% compute output & jacobian
        [a1, b1, tx, jac_a1_x, jac_a1_uv, jac_b1_x, jac_b1_uv, jac_tx_x, jac_tx_uv] = cons_a2(X, inde, localstamp, ta, tb, ind_1, h, param, Kinv, para);

        % output
        mm = b1'*tx*a1;
        if usehuber == 1
            [yhat(i), costd] = huber_cost(mm, delta);
        else
            yhat(i) = mm; costd = 1;
        end
        % jacobian
        cons1 = costd*(b1'*tx);
        cons2 = costd*(tx*a1)';
        cons3 = costd*(-b1'*skewm(a1));
        
        H(i,:) = cons1*jac_a1_x + cons2 * jac_b1_x + cons3 * jac_tx_x;
        
        Hf(i,:) = cons1*jac_a1_uv + cons2*jac_b1_uv + cons3 * jac_tx_uv;
        
        JRJt(i,i) = (cons1*jac_a1_uv)*(cons1*jac_a1_uv)' + (cons2*jac_b1_uv)*(cons2*jac_b1_uv)' + (cons3*jac_tx_uv)*(cons3*jac_tx_uv)';
        JRJt(i,i) = JRJt(i,i).*puv;
    end
    %% left null basis
%     [U,~,~] = svd(Hf);
%     AA = U(:,size(Hf,2)+1:end);
%     H = AA'*H;
%     JRJt = puv*eye(size(H,1));
%     yhat = AA'*yhat;
end

function [a, b, tx, jac_a_x, jac_a_uv, jac_b_x, jac_b_uv, jac_tx_x, jac_tx_uv] = cons_a2(X, inde, localstamp, ta, tb, ind_1, h, param, Kinv, para)
    fy = param.fy; fz = param.fz;
    fy1 = param.fy1; fz1 = param.fz1;
    dfy1 = param.dfy1; dfz1 = param.dfz1;
    nfy1 = param.nfy1; nfz1 = param.nfz1;
    r1 = param.r1; r2 = param.r2;
    distcorr1 = param.distcorr1; distcorr2 = param.distcorr2;
    k1 = param.k1;k2 = param.k2;k3 = param.k3;

    % call to get relative rotation
    [Ry, dty, iy, jy, R0y, R1y, xity, sy, Py, P0y, P1y] = geodesic_interpolate_rot_pos_vel(X, inde, localstamp, ta(ind_1));
    [Rz, dtz, iz, jz, R0z, R1z, xitz, sz, Pz, P0z, P1z] = geodesic_interpolate_rot_pos_vel(X, inde, localstamp, tb(ind_1));
    
    
    b = Rz*nfz1(:,ind_1);
    a = Ry*nfy1(:,ind_1);
    tx = skewm(Pz-Py);
%     tx = dtx ./ (norm(dtx)+1e-6); 
    
    xi1 = logSO3(R0y'*R1y);
    Jl1inv = leftJinv(xi1); Jr1inv = rightJinv(xi1);
    Jr2 = rightJ(xity);
    
    jac_R0 = Ry*skewm(nfy1(:,ind_1))*Jr2*Jl1inv*R0y'*dty-skewm(a);
    jac_R1 = -Ry*skewm(nfy1(:,ind_1))*Jr2*Jr1inv*R1y'*dty;
    
    tmp1 = -Ry*skewm(nfy1(:,ind_1))*xi1;
    
    %% jacobian
    jac_a_x = zeros(3,inde.cov_nongroup(end));
	if isfield(inde, 'cxy')
        ccx = para.cx*dft(X(inde.cxy(1)));
        ccy = para.cy*dft(X(inde.cxy(2)));
        ccf = para.f*dft(X(inde.f));
        
    	tmp2 = Ry*jac_f5(dfy1(:,ind_1));
        
        jac_cx = [-(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)^2/(fy1(3,ind_1)^2)-distcorr1(ind_1); ...
                  -(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*fy1(2,ind_1)/(fy1(3,ind_1)^2); ...
                  0];
        
        jac_cy = [-(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*fy1(2,ind_1)/(fy1(3,ind_1)^2); ...
                  -(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(2,ind_1)^2/(fy1(3,ind_1)^2)-distcorr1(ind_1); ...
                  0];
              
        jac_f = [-(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*(fy1(1,ind_1)/fy1(3,ind_1)^3+fy1(2,ind_1)/fy1(3,ind_1)^3); ...
                  -(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(2,ind_1)*(fy1(1,ind_1)/fy1(3,ind_1)^3+fy1(2,ind_1)/fy1(3,ind_1)^3); ...
                  1] + distcorr1(ind_1)*[0;0;1];
        
    	jac_a_x(:,[inde.cov_cxy, inde.cov_f]) = [tmp2*ccx*jac_cx tmp2*ccy*jac_cy tmp2*ccf*jac_f];
    
        if isfield(inde, 'k1')
            jac_a_x(:,[inde.cov_k1]) = tmp2*[r1(ind_1)*fy1(1,ind_1);r1(ind_1)*fy1(2,ind_1);0];
        end
        
        if isfield(inde, 'k2')
            jac_a_x(:,[inde.cov_k2]) = tmp2*[r1(ind_1)^2*fy1(1,ind_1);r1(ind_1)^2*fy1(2,ind_1);0];
        end
        
        if isfield(inde, 'k3')
            jac_a_x(:,[inde.cov_k3]) = tmp2*[r1(ind_1)^3*fy1(1,ind_1);r1(ind_1)^3*fy1(2,ind_1);0];
        end
    
        jac_u = [(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)^2/(fy1(3,ind_1)^2)+distcorr1(ind_1); ...
                  (2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*fy1(2,ind_1)/(fy1(3,ind_1)^2); ...
                  0];

        jac_v = [(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*fy1(2,ind_1)/(fy1(3,ind_1)^2); ...
                 (2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(2,ind_1)^2/(fy1(3,ind_1)^2)+distcorr1(ind_1); ...
                  0];
	else
    	tmp2 = Ry*jac_f5(dfy1(:,ind_1));
        if isfield(inde, 'k1')
            jac_a_x(:,[inde.cov_k1]) = tmp2*[r1(ind_1)*fy1(1,ind_1);r1(ind_1)*fy1(2,ind_1);0];
        end
        
        if isfield(inde, 'k2')
            jac_a_x(:,[inde.cov_k2]) = tmp2*[r1(ind_1)^2*fy1(1,ind_1);r1(ind_1)^2*fy1(2,ind_1);0];
        end
        
        if isfield(inde, 'k3')
            jac_a_x(:,[inde.cov_k3]) = tmp2*[r1(ind_1)^3*fy1(1,ind_1);r1(ind_1)^3*fy1(2,ind_1);0];
        end
        f1 = para.fx;
        f2 = para.fy;
        jac_u = [(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)^2/(f1)+distcorr1(ind_1)/f1; ...
                  (2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*fy1(2,ind_1)/(f1); ...
                  0];

        jac_v = [(2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(1,ind_1)*fy1(2,ind_1)/(f2); ...
                 (2*k1+4*k2*r1(ind_1)+6*k3*r1(ind_1)^2)*fy1(2,ind_1)^2/(f2)+distcorr1(ind_1)/f2; ...
                  0];
    end
    
    jac_a_uv(:,1:4) = [tmp2*jac_u tmp2*jac_v zeros(3,1) zeros(3,1)]+[zeros(3,1) tmp1*1/h*sy zeros(3,1) zeros(3,1)];

    v = fy(2,ind_1);
    ccc = para.ts*dft(X(inde.ts))*v/h;    
    ccd = para.td*dft(X(inde.td));    
    jac_a_x(:,[inde.cov_ts, inde.cov_td]) = [tmp1*sy*ccc tmp1*sy*ccd];
    
    jac_a_x(:,inde.cov_group(iy*9-8:iy*9-6)) = jac_R0;
    jac_a_x(:,inde.cov_group(jy*9-8:jy*9-6)) = jac_R1;
    
    %% b
    xi1 = logSO3(R0z'*R1z);
    Jl1inv = leftJinv(xi1); Jr1inv = rightJinv(xi1);
    Jr2 = rightJ(xitz);
    
    jac_R0 = Rz*skewm(nfz1(:,ind_1))*Jr2*Jl1inv*R0z'*dtz-skewm(b);
    jac_R1 = -Rz*skewm(nfz1(:,ind_1))*Jr2*Jr1inv*R1z'*dtz;
    
    tmp1 = -Rz*skewm(nfz1(:,ind_1))*xi1;
    
    jac_b_x = zeros(3,inde.cov_nongroup(end));
	if isfield(inde, 'cxy')
        ccx = para.cx*dft(X(inde.cxy(1)));
        ccy = para.cy*dft(X(inde.cxy(2)));
        ccf = para.f*dft(X(inde.f));
	    tmp2 = Rz*jac_f5(dfz1(:,ind_1));
        
        jac_cx = [-(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)^2/(fz1(3,ind_1)^2)-distcorr2(ind_1); ...
                  -(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*fz1(2,ind_1)/(fz1(3,ind_1)^2); ...
                  0];
        
        jac_cy = [-(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*fz1(2,ind_1)/(fz1(3,ind_1)^2); ...
                  -(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(2,ind_1)^2/(fz1(3,ind_1)^2)-distcorr2(ind_1); ...
                  0];
              
        jac_f = [-(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*(fz1(1,ind_1)/fz1(3,ind_1)^3+fz1(2,ind_1)/fz1(3,ind_1)^3); ...
                  -(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(2,ind_1)*(fz1(1,ind_1)/fz1(3,ind_1)^3+fz1(2,ind_1)/fz1(3,ind_1)^3); ...
                  1] + distcorr2(ind_1)*[0;0;1];
  
    	jac_b_x(:,[inde.cov_cxy, inde.cov_f]) = [tmp2*ccx*jac_cx tmp2*ccy*jac_cy tmp2*ccf*jac_f];
                                                   
        if isfield(inde, 'k1')
            jac_b_x(:,[inde.cov_k1]) = tmp2*[r2(ind_1)*fz1(1,ind_1);r2(ind_1)*fz1(2,ind_1);0];
        end
        
        if isfield(inde, 'k2')
            jac_b_x(:,[inde.cov_k2]) = tmp2*[r2(ind_1)^2*fz1(1,ind_1);r2(ind_1)^2*fz1(2,ind_1);0];
        end
        
        if isfield(inde, 'k3')
            jac_b_x(:,[inde.cov_k3]) = tmp2*[r2(ind_1)^3*fz1(1,ind_1);r2(ind_1)^3*fz1(2,ind_1);0];
        end
        
        jac_u = [(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)^2/(fz1(3,ind_1)^2)+distcorr2(ind_1); ...
                  (2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*fz1(2,ind_1)/(fz1(3,ind_1)^2); ...
                  0];

        jac_v = [(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*fz1(2,ind_1)/(fz1(3,ind_1)^2); ...
                 (2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(2,ind_1)^2/(fz1(3,ind_1)^2)+distcorr2(ind_1); ...
                  0];
    else
	    tmp2 = Rz*jac_f5(dfz1(:,ind_1));
        if isfield(inde, 'k1')
            jac_b_x(:,[inde.cov_k1]) = tmp2*[r2(ind_1);r2(ind_1);0];
        end
        
        if isfield(inde, 'k2')
            jac_b_x(:,[inde.cov_k2]) = tmp2*[r2(ind_1)^2;r2(ind_1)^2;0];
        end
        
        if isfield(inde, 'k3')
            jac_b_x(:,[inde.cov_k3]) = tmp2*[r2(ind_1)^3;r2(ind_1)^3;0];
        end
        
        f1 = para.fx;
        f2 = para.fy;
        jac_u = [(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)^2/(f1)+distcorr2(ind_1)/f1; ...
                  (2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*fz1(2,ind_1)/(f1); ...
                  0];

        jac_v = [(2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(1,ind_1)*fz1(2,ind_1)/(f2); ...
                 (2*k1+4*k2*r2(ind_1)+6*k3*r2(ind_1)^2)*fz1(2,ind_1)^2/(f2)+distcorr2(ind_1)/f2; ...
                  0];
    end
    
    jac_b_uv(:,1:4) = [ zeros(3,1) zeros(3,1) tmp2*jac_u tmp2*jac_v]+[zeros(3,1) zeros(3,1) zeros(3,1) tmp1*1/h*sz];
    
    v = fz(2,ind_1);
    ccc = para.ts*dft(X(inde.ts))*v/h;   
    ccd = para.td*dft(X(inde.td));  
    jac_b_x(:,[inde.cov_ts, inde.cov_td]) = [tmp1*sz*ccc tmp1*sz*ccd];
    
    jac_b_x(:,inde.cov_group(iz*9-8:iz*9-6)) = jac_R0;
    jac_b_x(:,inde.cov_group(jz*9-8:jz*9-6)) = jac_R1;
    
    jac_tx_x = zeros(3,inde.cov_nongroup(end));
    tmpz = (P1z - P0z);
    tmpy = (P1y - P0y);
    tr_z = para.ts*dft(X(inde.ts))*fz(2,ind_1)/h;   
    td_z = para.td*dft(X(inde.td));  
    
    tr_y = para.ts*dft(X(inde.ts))*fy(2,ind_1)/h;    
    td_y = para.td*dft(X(inde.td));    
    
    jac_tx_x(:,[inde.cov_ts, inde.cov_td]) = [tmpz*sz*tr_z-tmpy*sy*tr_y tmpz*sz*td_z-tmpy*sy*td_y];
    
    jac_tx_x(:,inde.cov_group(iz*9-5:iz*9-3)) = -dtz*eye(3);
    jac_tx_x(:,inde.cov_group(jz*9-5:jz*9-3)) = dtz*eye(3);
    jac_tx_x(:,inde.cov_group(iy*9-5:iy*9-3)) = dty*eye(3);
    jac_tx_x(:,inde.cov_group(jy*9-5:jy*9-3)) = -dty*eye(3);
    
    jac_tx_x(:,[inde.cov_ts, inde.cov_td]) = [tmpz*sz*tr_z-tmpy*sy*tr_y tmpz*sz*td_z-tmpy*sy*td_y];
    
    jac_tx_uv(:,1:4) = [zeros(3,1) -tmpy*sy*1/h zeros(3,1) tmpz*sz*1/h];
end

function varargout = geodesic_interpolate_rot_pos_vel(X, inde, localstamp, t)
% FUNCTION_NAME - geodesic interpolation of rotation matrices.
%
% Syntax:  varargout = geodesic_interpolate_rot(X, inde, localstamp, t)
%
% Inputs:
%    X - staets.
%    inde - index.
%    localstamp - imu timestamp.
%    t - feature timestamp.
%
% Outputs:
%    (1) - interpolated rotation.
%    (2) - relative time.
%    (3) - index of R0.
%    (4) - index of R1.
%    (5) - R0.
%    (6) - R1.
%    (7) - logSO3(R0'*R1)*dt;.
%    (8) - 1/(t_1-t_0).
%
% Author: Xiao Hu
% Technical University of Denmark
% email: xiahaa@space.dtu.dk
% Jan 2020; Last revision: 31-Jan-2020
%------------- BEGIN CODE --------------
    j = find(localstamp>t,1);
    i = j - 1;
    if i < 0, i = 1; warning('extra-interpolation!!! not accurate!!!'); end
    dt = (t - localstamp(i))/(localstamp(j)-localstamp(i));
    s = 1 / (localstamp(j)-localstamp(i));
    i = length(localstamp)+1 - i;
    j = length(localstamp)+1 - j;
    R0 = reshape(X(inde.group(i*15-14:i*15-6)),3,3);
    R1 = reshape(X(inde.group(j*15-14:j*15-6)),3,3);
    xit = logSO3(R0'*R1)*dt;
    Rt = R0*expSO3(xit);
    
    P0 = X(inde.group(i*15-5:i*15-3))';
    P1 = X(inde.group(j*15-5:j*15-3))';
    Pt = P0 + (P1-P0)*dt;
    
    varargout{1} = Rt;
    varargout{2} = dt;
    varargout{3} = i;
    varargout{4} = j;
    varargout{5} = R0;
    varargout{6} = R1;
    varargout{7} = xit;
    varargout{8} = s;
    varargout{9} = Pt;
    varargout{10} = P0;
    varargout{11} = P1;
end

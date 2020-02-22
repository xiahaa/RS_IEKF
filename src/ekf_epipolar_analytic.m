function [mean_est, var_est, npara] = ekf_epipolar_analytic(match_idx, match_x1, match_x2, gyrostamp, gyrogap, anglev, framestamp, para, endidx)
% This is the re-implementation of the EKF filter proposed in 
% Jia, Chao, and Brian L. Evans. 
% "Online camera-gyroscope autocalibration for cell phones." 
% IEEE Transactions on Image Processing 23.12 (2014): 5070-5081.
%
% The major reason for reimplement this method is that the original
% implementation is not available. 
% I reimplement this by carefully checking the method proposed in the paper.
% However, due to the limitation of my knowledge, there is no guarantee
% that my implementation does exactly the same as the proposed method.
% 
% Compare with the paper, there is an explicit change I made: I use the
% SO(3) rather than the unit-quaternion as the state for Ric. The reason
% for using SO(3) is that we do not need any reparameterization.

    addpath opt

    h = para.h;% img height
    w = para.w;% img width
    gyro_sigma = para.sigma;%
    pn = para.pn;% feature point noise
    maxtd = 0.03; % maximum range of parameter td is (+-) 0.03 s
    
    % initial variance of different parameters
    oc_sigma = 6.67^2;
    ts_sigma = 0.00167^2;
    td_sigma = 0.01^2;
    f_sigma = 20^2;
    bias_sigma = 0.0067^2;
    rcam_sigma = 0.0067^2;

    % initilize the EKF
    x_k_k = [para.cx; para.cy; para.ts; para.td; para.f; para.wd'; para.rcam];
    p_k_k = eye(11); 
    p_k_k(1,1) = oc_sigma; p_k_k(2,2) = oc_sigma; 
    p_k_k(3,3) = ts_sigma; p_k_k(4,4) = td_sigma; 
    p_k_k(5,5) = f_sigma; 
    p_k_k(6,6) = bias_sigma; p_k_k(7,7) = bias_sigma; p_k_k(8,8) = bias_sigma; 
    p_k_k(9,9) = rcam_sigma; p_k_k(10,10) = rcam_sigma; p_k_k(11,11) = rcam_sigma; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    match_start = 1;
    % function ekf_info_z finds relevant features from matched result (a txt
    % file) and return two vectors, matchstart indicates the corresponding line
    % number in the txt file.
    [z_k, match_start, y_k] = ekf_info_z (1, match_start, match_idx, match_x1, match_x2);
    % since each feature contains 2 position, so the feature number is the half
    % of the length of y_k
    fnum = length(y_k)/2;
    % find the gyro reading group
    idxstart = 1;
    % find_group function finds relevant gyro measurements' indices by using
    % the timestamp logged and delay (vision and gyro, should be small
    % enough). idxstart indicates the start searching id.
    [idxg, idxstart] = find_group(gyrostamp, framestamp(1) + para.td-maxtd , framestamp(3) + para.td+maxtd, idxstart);
    % roll back idxstart
    % this roll back is because gyro measurements may overlap with 2 frames.
    idxstart = idxstart - length(idxg) - 2;
    % this num is the augmented state used for angular velocity.
    x_num = length(idxg);

    % EKF predict
    x_k_km1 = x_k_k;
    p_k_km1 = p_k_k;   
    % EKF update
    if z_k ~= Inf
        % function epip_gen_hid generates random indices for epip constaints by
        % random perm 1:fnum and then reshape it as a Mx3 matrix.
        % hid = epip_gen_hid (fnum);% find hid randomly as proposed in
        % GlobalSip.
        hid = epip_gen_hid_analytical (y_k);% TIP

        [mm, H_k, VRV_k] = epip_info_meas_analytical (x_k_km1, hid, gyrostamp(idxg), gyrogap(idxg), ...
            framestamp(1), framestamp(2), z_k, y_k, h, reshape(anglev(idxg,:)',3*x_num,1), pn, gyro_sigma);

        % traditional KF update
        [x_k_k, p_k_k] = ekf_update_analytical (x_k_km1, p_k_km1, H_k, VRV_k, zeros(size(hid,1),1), mm);
    else % z_k == Inf means that there are no matched features detected
        x_k_k = x_k_km1;
        p_k_k = p_k_km1;
    end

    % log
    var_est(1,:) = [diag(p_k_k(1:11,1:11))'];
    mean_est(1,:) = x_k_k(1:11)';
    fprintf('idx = %d, ts = %f, td = %f, ocw = %f, och = %f, f = %f bias = %f %f %f rcamx = %f %f %f \n',1, x_k_k(3),x_k_k(4),x_k_k(1),x_k_k(2),x_k_k(5),x_k_k(6),x_k_k(7),x_k_k(8), x_k_k(9), x_k_k(10), x_k_k(11));
    % for the following every other pair of adjacent frames
    for i = 3:2:endidx
        tic;
        oldidxg = idxg;
        [idxg, idxstart] = find_group(gyrostamp, framestamp(i)+ para.td-maxtd, framestamp(i+2) + para.td+maxtd, idxstart);
        idxstart = idxstart - length(idxg) - 2;
        x_num = length(idxg);
        [z_k, match_start, y_k] = ekf_info_z (i, match_start, match_idx, match_x1, match_x2);
        fnum = length(y_k)/2;

        % EKF predict
        [x_k_km1, p_k_km1] = ekf_predict_epip_analytical (x_k_k, p_k_k);
        % EKF update
        if z_k ~= Inf
            % hid = epip_gen_hid (fnum);% GlobalSip
            hid = epip_gen_hid_analytical (y_k);% TIP
                    
            [mm, H_k, VRV_k] = epip_info_meas_analytical (x_k_km1, hid, gyrostamp(idxg), gyrogap(idxg), ...
                framestamp(i), framestamp(i+1), z_k, y_k, h, reshape(anglev(idxg,:)',3*x_num,1), pn, gyro_sigma);
            % traditional KF update
            [x_k_k, p_k_k] = ekf_update_analytical (x_k_km1, p_k_km1, H_k, VRV_k, zeros(size(hid,1),1), mm);
        else
            x_k_k = x_k_km1;
            p_k_k = p_k_km1;
        end
        % make sure state variables satisfy constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if x_k_k(3) < 0
            x_k_k(3) = 0;
        elseif x_k_k(3) >(1/30)
            x_k_k(3) = 1/30;
        end
        if x_k_k(4) < para.td-maxtd
            x_k_k(4) = para.td-maxtd;
        elseif x_k_k(4) > para.td+maxtd
            x_k_k(4) = para.td+maxtd;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % log
        var_est(i,:) = [diag(p_k_k(1:11,1:11))'];
        mean_est(i,:) = x_k_k(1:11)';
        fprintf('idx = %d, ts = %f, td = %f, ocw = %f, och = %f, f = %f bias = %f %f %f rcamx = %f %f %f \n',i, x_k_k(3),x_k_k(4),x_k_k(1),x_k_k(2),x_k_k(5),x_k_k(6),x_k_k(7),x_k_k(8), x_k_k(9), x_k_k(10), x_k_k(11));
        runtime(i) = toc;

    end
    npara.ts = x_k_k(3);
    npara.td = x_k_k(4);
    npara.wd = [x_k_k(6) x_k_k(7) x_k_k(8)];
    npara.f = x_k_k(5);
    npara.cx = x_k_k(1);
    npara.cy = x_k_k(2);
    npara.rcam = [x_k_k(9); x_k_k(10); x_k_k(11)];

    disp(['mean runtime is ', num2str(mean(runtime(runtime ~= 0)))]);

end

function hid = epip_gen_hid_analytical (y_k)
    % generate groups of features (3 in each group) for the coplanar
    % constraint according to the principle described in the TIP paper.
    % Jia, Chao, and Brian L. Evans. 
    % "Online camera-gyroscope autocalibration for cell phones." 
    % IEEE Transactions on Image Processing 23.12 (2014): 5070-5081.
    y_k = reshape(y_k, 2, length(y_k)/2);
    [~,inds] = sort(y_k(2,:));
    fnum = floor(size(y_k,2)/3)*3;
    hid = reshape(inds(1:fnum),fnum/3,3);    
end

function [x_k_km1, p_k_km1] = ekf_predict_epip_analytical (x_k_k, p_k_k)
    % ekf filtering prediction step
    % start prediction
    x_k_km1 = x_k_k;
    p_k_km1 = p_k_k;
end

function [x_k_k, p_k_k] = ekf_update_analytical (x_k_km1, p_k_km1, H_k, VRV_k, z_k, py)
    % EKF filtering update step
    K = p_k_km1*H_k'/(H_k*p_k_km1*H_k'+VRV_k); % ekf gain matrix
    x_k_k = [x_k_km1(1:end-3);x_k_km1(end-2:end)] + K*(z_k - py(:));
    p_k_k = (eye(length(x_k_k)) - K*H_k)*p_k_km1;
end
	
function [yhat, H, JRJt] = epip_info_meas_analytical(x, hid, localstamp, localgap, fta, ftb, fz, fy, h, anglev, puv, pw)
    % split current states
    ocx = x(1);
    ocy = x(2);
    tr = x(3);
    td = x(4);
    f = x(5);
    bias = x(6:8);
    rcam = x(9:11);

    Q = epip_get_Q(rcam);% this Q is the relative rotation from imu to camera

    hnum = size(hid,1);
    H = zeros(hnum, length(x));
    yhat = zeros(hnum,1);
    
    fta = fta + td;      % aligned frame time
    ftb = ftb + td; 
    ta = fta + fy(2:2:end)./h.*tr;% aligned feature time
    tb = ftb + fz(2:2:end)./h.*tr;
    tref = ftb;          % epipolar constraint establishing time
    
    % since for rs camera, each feature will have its own rotation matrix
    % because they may lay on different row. z means features in new frame, y
    % means features in old frame. 
    
    %% to d
    fy1 = [reshape(fy, 2, length(fy)/2); f*ones(1, length(fy)/2)]; fy1 = fy1 - [ocx;ocy;0];
    fz1 = [reshape(fz, 2, length(fz)/2); f*ones(1, length(fz)/2)]; fz1 = fz1 - [ocx;ocy;0];
    
    %% normalization to e
    nfy1 = fy1 ./ vecnorm(fy1);
    nfz1 = fz1 ./ vecnorm(fz1);
    
    %% to imu
    nfy1_imu = Q' * nfy1;
    nfz1_imu = Q' * nfz1;
    
    %% angular velocities after bias correction.
    w_num = length(localstamp);
    anglev = anglev + repmat(bias, w_num, 1);
    
    param = struct;
    param.fy = reshape(fy, 2, length(fy)/2); param.fz = reshape(fz, 2, length(fz)/2);
    param.fy1 = fy1; param.fz1 = fz1;
    param.nfy1 = nfy1; param.nfz1 = nfz1;
    param.nfy1_imu = nfy1_imu; param.nfz1_imu = nfz1_imu;
    
    JRJt1 = zeros(hnum, hnum);
    JRJt2 = zeros(hnum, hnum);
    
    for i = 1:hnum
        %% compute intermediate variables 
        % get id, the combinations (groups of 3 features) are in hid
        ind_1 = hid(i,1);
        ind_2 = hid(i,2);
        ind_3 = hid(i,3);
        
        %% compute output & jacobian
        [a1, b1, jac_a1_x, jac_a1_uv, jac_a1_w, jac_b1_x, jac_b1_uv, jac_b1_w] = cons_a(x, localstamp, localgap, ...
            ta, tb, tref, w_num, ind_1, anglev, h, param, Q);

        [a2, b2, jac_a2_x, jac_a2_uv, jac_a2_w, jac_b2_x, jac_b2_uv, jac_b2_w] = cons_a(x, localstamp, localgap, ...
            ta, tb, tref, w_num, ind_2, anglev, h, param, Q);        
        
        [a3, b3, jac_a3_x, jac_a3_uv, jac_a3_w, jac_b3_x, jac_b3_uv, jac_b3_w] = cons_a(x, localstamp, localgap, ...
            ta, tb, tref, w_num, ind_3, anglev, h, param, Q);        

        % output
        f1 = cross(a1,b1); f2 = cross(a2,b2); f3 = cross(a3,b3); 
        yhat(i) = dot(cross(f1,f2),f3);
        % jacobian
        cons1 = jac_f3(f1,f2,f3,1);
        cons2 = jac_f3(f1,f2,f3,2);
        cons3 = jac_f3(f1,f2,f3,3);
        
        tmp1 = cons1 * jac_f4(a1,b1,1);
        tmp2 = cons1 * jac_f4(a1,b1,2);
        tmp3 = cons2 * jac_f4(a2,b2,1);
        tmp4 = cons2 * jac_f4(a2,b2,2);
        tmp5 = cons3 * jac_f4(a3,b3,1);
        tmp6 = cons3 * jac_f4(a3,b3,2);
        
        H(i,:) = tmp1*jac_a1_x + tmp2*jac_b1_x + ...
                 tmp3*jac_a2_x + tmp4*jac_b2_x + ...
                 tmp5*jac_a3_x + tmp6*jac_b3_x;
        
        JRJt1(i,i) = (tmp1*jac_a1_uv*(tmp1*jac_a1_uv)' + tmp2*jac_b1_uv*(tmp2*jac_b1_uv)' +...
                      tmp3*jac_a2_uv*(tmp3*jac_a2_uv)' + tmp4*jac_b2_uv*(tmp4*jac_b2_uv)' +...
                      tmp5*jac_a3_uv*(tmp5*jac_a3_uv)' + tmp6*jac_b3_uv*(tmp6*jac_b3_uv)').*puv;
                 
        JRJt2(i,i) = (tmp1*jac_a1_w*(tmp1*jac_a1_w)' + tmp2*jac_b1_w*(tmp2*jac_b1_w)' +...
                      tmp3*jac_a2_w*(tmp3*jac_a2_w)' + tmp4*jac_b2_w*(tmp4*jac_b2_w)' +...
                      tmp5*jac_a3_w*(tmp5*jac_a3_w)' + tmp6*jac_b3_w*(tmp6*jac_b3_w)').*pw;
    end
    JRJt = JRJt1 + JRJt2;
end

function [a, b, jac_a_x, jac_a_uv, jac_a_w, jac_b_x, jac_b_uv, jac_b_w] = cons_a(x, localstamp, localgap, ta, tb, tref, w_num, ind_1, anglev, h, param, Q)
    fy = param.fy; fz = param.fz;
    fy1 = param.fy1; fz1 = param.fz1;
    nfy1_imu = param.nfy1_imu; nfz1_imu = param.nfz1_imu;

    % call to get relative rotation
    dty = find_dt (localstamp, localgap, ta(ind_1), tref, w_num);
    [Ry, Rys, nzy, wys] = epip_R_gen_analytical (anglev, dty, 1);% backward
    
    dtz = find_dt (localstamp, localgap, tref, tb(ind_1), w_num);
    [Rz, Rzs, nzz, wzs] = epip_R_gen_analytical (anglev, dtz, 2);% forward

    %
    a = Q*Ry*nfy1_imu(:,ind_1);
    b = Q*Rz*nfz1_imu(:,ind_1);
    
    %% jacobian
    jac_a_x = zeros(3,length(x));
    tmp1 = Q*Ry*Q'*jac_f5(fy1(:,ind_1));
    jac_a_x(:,[1,2,5]) = [tmp1*[-1;0;0] tmp1*[0;-1;0] tmp1*[0;0;1]];
    jac_a_uv(:,1:2) = [tmp1*[1;0;0] tmp1*[0;1;0]];
    
    [tmp2,tmp3] = jac_f8(Rys, wys, nfy1_imu(:,ind_1), fy(:,ind_1), h, 1);
    tmp2 = Q * tmp2; tmp3 = Q * tmp3;
    jac_a_x(:,[4,3]) = [tmp2 tmp3];
    
    [tmp4, tmp5] = jac_f9(Rys, dty(nzy), nfy1_imu(:,ind_1),1);
    tmp4 = Q * tmp4;
    jac_a_x(:,6:8) = tmp4;
    
    jac_a_w = zeros(3, 3*length(dty));
    for i = 1:length(nzy)
        jac_a_w(:,3*nzy(i)-2:3*nzy(i)) = Q * tmp5(:,:,i);
    end
    
    Jl = rightJ(logSO3(Q));
    jac_a_x(:,9:11) = -Q*skewm(Ry*nfy1_imu(:,ind_1))*Jl + Q*Ry*skewm(nfy1_imu(:,ind_1))*Jl;
    
    %% b
    jac_b_x = zeros(3,length(x));
    tmp1 = Q*Rz*Q'*jac_f5(fz1(:,ind_1));
    jac_b_x(:,[1,2,5]) = [tmp1*[-1;0;0] tmp1*[0;-1;0] tmp1*[0;0;1]];
    jac_b_uv(:,1:2) = [tmp1*[1;0;0] tmp1*[0;1;0]];
    
    [tmp2,tmp3] = jac_f8(Rzs, wzs, nfz1_imu(:,ind_1), fz(:,ind_1), h, 2);
    tmp2 = Q * tmp2; tmp3 = Q * tmp3;% no minus because of forward 
    jac_b_x(:,[4,3]) = [tmp2 tmp3];
    
    [tmp4, tmp5] = jac_f9(Rzs, dtz(nzz), nfz1_imu(:,ind_1),2);
    tmp4 = Q * tmp4;
    jac_b_x(:,6:8) = tmp4;
    
    jac_b_w = zeros(3, 3*length(dtz));
    for i = 1:length(nzz)
        jac_b_w(:,3*nzz(i)-2:3*nzz(i)) = Q * tmp5(:,:,i);
    end
    jac_b_x(:,9:11) = -Q*skewm(Rz*nfz1_imu(:,ind_1))*Jl + Q*Rz*skewm(nfz1_imu(:,ind_1))*Jl;
end

function varargout = epip_R_gen_analytical (angles, dt, opt)
    % compute the relative rotation matrix correspond to a certain feature
    % point
    a_num = length(angles)/3;
    angles = reshape(angles, 3, a_num)';

    nz = find(dt~=0);
    R = eye(3);
    Rs = zeros(3,3,a_num);
    
    if opt == 1
        % backward, actually is exp(-w*dt)
        for i = nz(end):-1:nz(1)
            omega = angles(i,:).*(-dt(i));
            theta = norm(omega);
            omega = omega./theta;
            skew_matrix = skewm(omega);
            Rlocal = eye(3) + sin(theta).*skew_matrix + (1-cos(theta)).*(skew_matrix*skew_matrix);
            R = R*Rlocal;
            Rs(:,:,i) = Rlocal;
        end
    else
        % forward, exp(w*dt)
        for i = nz(1):1:nz(end)
            omega = angles(i,:).*(dt(i));
            theta = norm(omega);
            omega = omega./theta;
            skew_matrix = skewm(omega);
            Rlocal = eye(3) + sin(theta).*skew_matrix + (1-cos(theta)).*(skew_matrix*skew_matrix);
            R = R*Rlocal;
            Rs(:,:,i) = Rlocal;
        end
    end
    varargout{1} = R;% final R
    varargout{2} = Rs(:,:,nz);% used for jacobian
    varargout{3} = nz;% used for jacobian
    varargout{4} = angles(nz,:)';% used for jacobian
end

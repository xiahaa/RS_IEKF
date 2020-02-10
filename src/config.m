function para = config()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter type:
%   1. EKF proposed in (use author-provided source code)
%       Jia, Chao, and Brian L. Evans. 
%       "Online calibration and synchronization of cellphone camera and gyroscope." 
%       2013 IEEE Global Conference on Signal and Information Processing. IEEE, 2013.
%   2. EKF proposed in (re-implementation)
%       Jia, Chao, and Brian L. Evans. 
%       "Online camera-gyroscope autocalibration for cell phones." 
%       IEEE Transactions on Image Processing 23.12 (2014): 5070-5081.
%   3. IEKF 
%
    flt_type = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if flt_type == 1 || flt_type == 2
        % set the initial values of parameters
        para.ts = 0.025; % readout time of rolling shutter camera, 0.025 
        para.td = 1.175;%1.175; % timestamp delay	
        para.wd = [0 0 0]; % gyroscope bias
        para.f = 700; % camera focal length
        para.w = 720; % frame width (this is a fixed value)
        para.h = 480; % frame height (this is a fixed value)
        para.cx = (para.w-1)/2; % principal point
        para.cy = (para.h-1)/2; % principal point
        para.rcam = [pi/sqrt(2); -pi/sqrt(2); 0]; % relative orientation between camera and gyro

        para.sigma = 1e-6; % the noise variance of gyro readings
        para.pn = 1; % noise variance of feature detection
    else
        % set the initial values of parameters
        para.ts = 0.025; % readout time of rolling shutter camera, 0.025 
        para.td = 1.175;%1.175; % timestamp delay	
        para.wd = [0 0 0]; % gyroscope bias
        para.f = 700; % camera focal length
        para.fx = para.f; para.fy = para.f; 
        para.w = 720; % frame width (this is a fixed value)
        para.h = 480; % frame height (this is a fixed value)
        para.cx = (para.w-1)/2; % principal point
        para.cy = (para.h-1)/2; % principal point
        para.rcam = [pi/sqrt(2); -pi/sqrt(2); 0]; % relative orientation between camera and gyro

        para.sigma = 1e-6; % the noise variance of gyro readings
        para.pn = 1; % noise variance of feature detection
    end
    
    if flt_type == 1
        para.flt = @ekf_epipolar;
    elseif flt_type == 2
        para.flt = @ekf_epipolar_analytic;
    elseif flt_type == 3
        para.flt = @InEKF_rs_calib_rep;
    else
        error('Wrong filter type!!!!');
    end
end
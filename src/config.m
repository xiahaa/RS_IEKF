function para = config(varargin)
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
    if flt_type == 1
        para.flt = @ekf_epipolar;
    elseif flt_type == 2
        para.flt = @ekf_epipolar_analytic;
    elseif flt_type == 3
        para.flt = @InEKF_rs_calib_rep;
    else
        error('Wrong filter type!!!!');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isempty(varargin) || varargin{1} == 1
        % use Jia Chao's dataset
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
    else
        paths = {'../data/rccar',...
                 '../data/rotation',...
                 '../data/walk'};
        id = 2;
        td = [5.2, 2.8, 3.39];% GT = [5.208, 2.838, 3.3858]; init around GT with some errors.

        para.td = td(id);%;%2.83;%1.175; % timestamp delay	
        % set the initial values of parameters
        para.ts = 0.033; % readout time of rolling shutter camera
        para.wd = [0 0 0]; % gyroscope bias
        para.fx = 853.12703455; % camera focal length
        para.fy = 873.54956631; % camera focal length
        para.f = (para.fx+para.fy)/2; % camera focal length
        % para.fx = 863.3383; % camera focal length
        para.w = 1920; % frame width (this is a fixed value)
        para.h = 1080; % frame height (this is a fixed value)
        para.cx = 988.06311256; % principal point
        para.cy = 525.71056312; % principal point
        para.rcam = [0;0;-1.54]; % relative orientation between camera and gyro

        para.sigma = 5*1e-3; % the noise variance of gyro readings
        para.pn = 1; % noise variance of feature detection
        para.basepath = paths{id};
    end    
end
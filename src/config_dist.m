function para = config_dist(varargin)
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
        para.flt = @InEKF_rs_calib_rep_exp;
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
            para.dist = [0,0,0];% k1 k2 k3 only radial distortion
            
            para.sigma = 1e-6; % the noise variance of gyro readings
            para.pn = 1; % noise variance of feature detection
        end
    else
        paths = {'../data/rccar',...
                 '../data/rotation',...
                 '../data/walk'};
             
        id = 1;
        td = [5.2, 2.8, 3.39];% GT = [5.208, 2.838, 3.3858]; init around GT with some errors.
        
        ric_truth = {[-0.0712 0.9974 -0.0130; -0.9954 -0.0719 -0.0640; -0.0648 0.0084 0.9979], ...
                     [ 0.0288 0.9996 -0.0011; -0.9996 0.0288 -0.0030;-0.0029 0.0012 1.0000], ...
                     [ 0.0213 0.9998 0.0024; -0.9996 0.0212 0.0164; 0.0164 -0.0028 0.9999]};
        
        tr_truth = [ 0.0316734,  0.0316734,  0.0316734];
        td_truth = [5.20867366868152E+00, 2.83818847389866E+00, 3.38581144830956E+00];
                 
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
        para.rcam = [0;0;-pi/2]; % relative orientation between camera and gyro
        para.dist = [0,0,0];% k1 k2 k3 only radial distortion
        para.sigma = 5*1e-3; % the noise variance of gyro readings
        para.pn = 1; % noise variance of feature detection
        para.basepath = paths{id};        
        para.fix = 1;
        para.ric_th = ric_truth{id};
        para.tr_th = tr_truth(id);
        para.td_th = td_truth(id);
    end    
end
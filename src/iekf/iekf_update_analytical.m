function [X, P] = iekf_update_analytical (X, inde, P, H_k, VRV_k, z_k, py, varargin)
%FUNCTION_NAME - iekf update.
%
% Syntax:  [X, P] = iekf_update_analytical (X, inde, P, H_k, VRV_k, z_k, py)
%
% Inputs:
%    X - states
%    inde - index
%    P - covariance
%    H_k - measurement Jacobian matrix.
%    VRV_k - measurement Cov matrix (implicit EKF).
%    z_k - real measurement.
%    py - predicted measurement.
%
% Outputs:
%    X - states
%    P - covariance
%
% Author: Xiao Hu
% Technical University of Denmark
% email: xiahaa@space.dtu.dk
% Jan 2020; Last revision: 31-Jan-2020
%------------- BEGIN CODE --------------
    add_pos_vel = 0;
    if ~isempty(varargin)
        add_pos_vel = varargin{1};
    end
    
    % EKF filtering update step
    S = H_k*P*H_k'+VRV_k;
    K = P*H_k'/(S); % ekf gain matrix
    residual = z_k - py(:);
    innovation = K*residual;
    
    %% group update
    if add_pos_vel == 0
        for i = 1:inde.group_num
            R = reshape(X(inde.group(i*9-8:i*9)),3,3);
            R = expSO3(innovation(inde.cov_group(i*3-2:i*3))) * R;
            X(inde.group(i*9-8:i*9)) = (R(:))';
        end
    else
        for i = 1:inde.group_num
            R = reshape(X(inde.group(i*15-14:i*15-6)),3,3);
            R = expSO3(innovation(inde.cov_group(i*9-8:i*9-6))) * R;
            X(inde.group(i*15-14:i*15-6)) = (R(:))';
            X(inde.group(i*15-5:i*15-3)) = X(inde.group(i*15-5:i*15-3)) + innovation(inde.cov_group(i*9-5:i*9-3))';
            X(inde.group(i*15-2:i*15)) = X(inde.group(i*15-2:i*15-0)) + innovation(inde.cov_group(i*9-2:i*9-0))';
        end
    end
	if isfield(inde, 'cxy')
    	X([inde.cxy, inde.f, inde.ts, inde.td, inde.wd]) = X([inde.cxy, inde.f, inde.ts, inde.td, inde.wd]) + innovation([inde.cov_cxy, inde.cov_f, inde.cov_ts, inde.cov_td, inde.cov_wd])';
	else
		X([inde.ts, inde.td, inde.wd]) = X([inde.ts, inde.td, inde.wd]) + innovation([inde.cov_ts, inde.cov_td, inde.cov_wd])';
    end
    
    if isfield(inde, 'k1')
        X([inde.k1]) = X([inde.k1]) + innovation([inde.cov_k1])';
    end
    if isfield(inde, 'k2')
        X([inde.k2]) = X([inde.k2]) + innovation([inde.cov_k2])';
    end
    if isfield(inde, 'k3')
        X([inde.k3]) = X([inde.k3]) + innovation([inde.cov_k3])';
    end
    
    R = reshape(X(inde.rcam),3,3);
    R = expSO3(innovation(inde.cov_rcam)) * R;
    X(inde.rcam) = (R(:))';
    P = P - K*S*K';%(eye(size(P),1) - K*H_k)*P;
end
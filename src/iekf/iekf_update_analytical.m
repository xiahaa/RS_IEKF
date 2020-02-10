function [X, P] = iekf_update_analytical (X, inde, P, H_k, VRV_k, z_k, py)
    % EKF filtering update step
    S = H_k*P*H_k'+VRV_k;
    K = P*H_k'/(S); % ekf gain matrix
    residual = z_k - py(:);
    innovation = K*residual;
    
    %% group update
    for i = 1:inde.group_num
        R = reshape(X(inde.group(i*9-8:i*9)),3,3);
        R = expSO3(innovation(inde.cov_group(i*3-2:i*3))) * R;
        X(inde.group(i*9-8:i*9)) = (R(:))';
    end
	if isfield(inde, 'cxy')
    	X([inde.cxy, inde.f, inde.ts, inde.td, inde.wd]) = X([inde.cxy, inde.f, inde.ts, inde.td, inde.wd]) + innovation([inde.cov_cxy, inde.cov_f, inde.cov_ts, inde.cov_td, inde.cov_wd])';
	else
		X([inde.ts, inde.td, inde.wd]) = X([inde.ts, inde.td, inde.wd]) + innovation([inde.cov_ts, inde.cov_td, inde.cov_wd])';
	end
    R = reshape(X(inde.rcam),3,3);
    R = expSO3(innovation(inde.cov_rcam)) * R;
    X(inde.rcam) = (R(:))';
    P = P - K*S*K';%(eye(size(P),1) - K*H_k)*P;
end
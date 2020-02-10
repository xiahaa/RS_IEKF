function [yhat, H, JRJt] = rep_info_meas_analytical(X, inde, localstamp, fta, ftb, fz, fy, h, puv,param)
%FUNCTION_NAME - analytical solution for measurement update.
%
% Syntax:  [yhat, H, JRJt] = rep_info_meas_analytical(X, inde, localstamp, fta, ftb, fz, fy, h, puv,param)
%
% Inputs:
%    X - states
%    inde - index
%    localstamp - imu timestamp
%    fta - timstamp of u
%    ftb - timestamp of u'
%    fz - u'
%    fy - u
%    h - img height
%    puv - cov of feature
%    param - auxiliary parameters.
%
% Outputs:
%    yhat - predicted measurement.
%    H - measurement Jacobian matrix.
%    JRJt - measurement Cov matrix (implicit EKF).
%
%
% Author: Xiao Hu
% Technical University of Denmark
% email: xiahaa@space.dtu.dk
% Jan 2020; Last revision: 31-Jan-2020
%------------- BEGIN CODE --------------
    x_nongroup = X(inde.nongroup);
    % split current states
	if isfield(inde,'cxy')
	    ocx = x_nongroup(1);
	    ocy = x_nongroup(2);
	    tr = x_nongroup(3);
	    td = x_nongroup(4);
	    f = x_nongroup(5);
	else
		tr = x_nongroup(1);
    	td = x_nongroup(2);
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
		K = [param.fx 0 param.cx;0 param.fy param.cy;0 0 1];Kinv = inv(K);
    	fy1 = [reshape(fy, 2, length(fy)/2); f*ones(1, length(fy)/2)]; fy1 = fy1 - [ocx;ocy;0];
    	fz1 = [reshape(fz, 2, length(fz)/2); f*ones(1, length(fz)/2)]; fz1 = fz1 - [ocx;ocy;0];
    else
		K = [param.fx 0 param.cx;0 param.fy param.cy;0 0 1];Kinv = inv(K);
    	fy1 = [reshape(fy, 2, length(fy)/2); ones(1, length(fy)/2)]; fy1 = Kinv * fy1;
    	fz1 = [reshape(fz, 2, length(fz)/2); ones(1, length(fz)/2)]; fz1 = Kinv * fz1;
	end
    %% normalization to e
    nfy1 = fy1 ./ vecnorm(fy1);
    nfz1 = fz1 ./ vecnorm(fz1);
    
    param = struct;
    param.fy = reshape(fy, 2, length(fy)/2); param.fz = reshape(fz, 2, length(fz)/2);
    param.fy1 = fy1; param.fz1 = fz1;
    param.nfy1 = nfy1; param.nfz1 = nfz1;
    
	Hf = zeros(hnum, 2);
    JRJt = zeros(hnum, hnum);
    delta = 0.05;
    usehuber = 1;
    for i = 1:hnum
        %% compute intermediate variables 
        ind_1 = i;
        
        %% compute output & jacobian
        [a1, b1, jac_a1_x, jac_a1_uv, jac_b1_x, jac_b1_uv] = cons_a(X, inde, localstamp, ta, tb, ind_1, h, param, Kinv);

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
end

function [a, b, jac_a_x, jac_a_uv, jac_b_x, jac_b_uv] = cons_a(X, inde, localstamp, ta, tb, ind_1, h, param, Kinv)
    fy = param.fy; fz = param.fz;
    fy1 = param.fy1; fz1 = param.fz1;
    nfy1 = param.nfy1; nfz1 = param.nfz1;

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
    	tmp2 = Ry*jac_f5(fy1(:,ind_1));
    	jac_a_x(:,[inde.cov_cxy, inde.cov_f]) = [tmp2*[-1;0;0] tmp2*[0;-1;0] tmp2*[0;0;1]];
	else
		tmp2 = Ry*jac_f5(fy1(:,ind_1))*diag([Kinv(1,1),Kinv(2,2),0]);
	end
    jac_a_uv(:,1:2) = [tmp2*[1;0;0] tmp2*[0;1;0]]+[zeros(3,1) tmp1*1/h*sy];

    jac_a_x(:,[inde.cov_ts, inde.cov_td]) = [tmp1*fy(2,ind_1)/h*sy tmp1*sy];
    
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
	    tmp2 = Rz*jac_f5(fz1(:,ind_1));
    	jac_b_x(:,[inde.cov_cxy, inde.cov_f]) = [tmp2*[-1;0;0] tmp2*[0;-1;0] tmp2*[0;0;1]];
    else
		tmp2 = Rz*jac_f5(fz1(:,ind_1))*diag([Kinv(1,1),Kinv(2,2),0]);
	end
	jac_b_uv(:,1:2) = [tmp2*[1;0;0] tmp2*[0;1;0]] + [zeros(3,1) tmp1*1/h*sz];

    jac_b_x(:,[inde.cov_ts, inde.cov_td]) = [tmp1*fz(2,ind_1)/h*sz tmp1*sz];
    
    jac_b_x(:,inde.cov_group(iz*3-2:iz*3)) = jac_R0;
    jac_b_x(:,inde.cov_group(jz*3-2:jz*3)) = jac_R1;
end
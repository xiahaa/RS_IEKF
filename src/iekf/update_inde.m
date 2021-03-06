function inde = update_inde(inde, x_k_k, p_k_k)
%FUNCTION_NAME - update index.
%
% Syntax:  inde = update_inde(inde, x_k_k, p_k_k)
%
% Inputs:
%    inde - old index
%    x_k_k - state
%    p_k_k - covariance
%
% Outputs:
%    inde - new index
%
%
% Author: Xiao Hu
% Technical University of Denmark
% email: xiahaa@space.dtu.dk
% Jan 2020; Last revision: 31-Jan-2020
%------------- BEGIN CODE --------------
    inde.group = 1:9*inde.group_num;
    inde.cov_group = 1:3*inde.group_num;
    
    inde.nongroup = (inde.group(end)+1):(inde.group(end)+length(x_k_k));
	if length(x_k_k) == 17
        % free cx, cy, f
    	inde.cxy = inde.nongroup(1:2);inde.f = inde.nongroup(5);
	    inde.ts = inde.nongroup(3);inde.td = inde.nongroup(4);
	    inde.wd = inde.nongroup(6:8);inde.rcam = inde.nongroup(9:17);
    
	    inde.cov_nongroup = (inde.cov_group(end)+1):(inde.cov_group(end)+size(p_k_k,1));
	    inde.cov_cxy = inde.cov_nongroup(1:2);inde.cov_f = inde.cov_nongroup(5);
	    inde.cov_ts = inde.cov_nongroup(3);inde.cov_td = inde.cov_nongroup(4);
	    inde.cov_wd = inde.cov_nongroup(6:8);inde.cov_rcam = inde.cov_nongroup(9:11);
    else
        % fix cx, cy, f
		inde.ts = inde.nongroup(1);inde.td = inde.nongroup(2);
    	inde.wd = inde.nongroup(3:5);inde.rcam = inde.nongroup(6:14);
    
    	inde.cov_nongroup = (inde.cov_group(end)+1):(inde.cov_group(end)+size(p_k_k,1));
    	inde.cov_ts = inde.cov_nongroup(1);inde.cov_td = inde.cov_nongroup(2);
    	inde.cov_wd = inde.cov_nongroup(3:5);inde.cov_rcam = inde.cov_nongroup(6:8);
	end
end
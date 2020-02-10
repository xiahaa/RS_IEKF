function [jac_td, jac_tr] = jac_f8(rs, ws, fi, z, h, opt)
    % jacobian of IIrs(tr,td): implicit function of tr, td
    %% for i = 1
    N = size(rs,3);
    R1 = eye(3);
    
    if opt == 2
        for i = 2:N
            R1 = R1*rs(:,:,i);
        end
        fi1 = R1 * fi;
        jac1 = -rs(:,:,1) * skewm(fi1);
    
        %% for i = M
        jac2 = -rs(:,:,1) * R1 * skewm(fi);
        
        jac_td = -jac1 * ws(:,1) + jac2 * ws(:, end);
        jac_tr = jac2 * ws(:, end) * (z(2)/h);
    else
        for i = N-1:-1:1
            R1 = R1*rs(:,:,i);
        end
        fi1 = R1 * fi;
        jac1 = -rs(:,:,N) * skewm(fi1);
    
        %% for i = M
        jac2 = -rs(:,:,N) * R1 * skewm(fi);
        
        jac_td = -jac1 * ws(:,end) + jac2 * ws(:, 1);
        jac_tr = jac2 * ws(:, 1) * (z(2)/h);
    end
end
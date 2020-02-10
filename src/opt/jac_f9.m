function [jac, jacs] = jac_f9(rs, ts, fi, opt)
    % jacobian of IIrs(bg): implicit function of bg
    jac = zeros(3,3);
    N = size(rs,3);
    R1 = eye(3);
    
    if opt == 2
        %% forward
        for i = 1:N
            R1 = R1 * rs(:,:,i);
        end
        f1 = fi;
        jacs = [];
        if nargout >= 2
            jacs = zeros(3,3,N);
            for  i = N:-1:1
                jacs(:,:,i) = (-ts(i).*(R1*skewm(f1)));
                jac = jac + jacs(:,:,i);
                R1 = R1 * rs(:,:,i)';
                f1 = rs(:,:,i)*f1;
            end
        else
            for  i = N:-1:1
                jac = jac + (-ts(i).*(R1*skewm(f1)));
                R1 = R1 * rs(:,:,i)';
                f1 = rs(:,:,i)*f1;
            end
        end
    else
        %% backward
        for i = N:-1:1
            R1 = R1*rs(:,:,i);
        end
        f1 = fi;
        jacs = [];
        if nargout >= 2
            jacs = zeros(3,3,N);
            for  i = 1:N
                jacs(:,:,i) = (ts(i).*(R1*skewm(f1)));
                jac = jac + jacs(:,:,i);
                R1 = R1 * rs(:,:,i)';
                f1 = rs(:,:,i)*f1;
            end
        else
            for  i = 1:N
                jac = jac + (ts(i).*(R1*skewm(f1)));
                R1 = R1 * rs(:,:,i)';
                f1 = rs(:,:,i)*f1;
            end
        end
    end
end
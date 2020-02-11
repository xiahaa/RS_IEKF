
function r = logSO3(R)
    r = [0;0;0];
    res = (trace(R)-1)*0.5;
    if res < -1
        res = -1;
    end
    angle = acos(res);
    if angle > 1e-10
        so3 = angle / (2*sin(angle))*(R-R');
        r = [-so3(2,3);so3(1,3);-so3(1,2)];
    else
        so3 = zeros(3,3);
        for i = 1:2
            so3 = so3 + (-1)^(i-1)/i.*(R-eye(3))^(i);
        end
        so3 = 0.5.*(so3-so3');
        r = [-so3(2,3);so3(1,3);-so3(1,2)];
    end
end
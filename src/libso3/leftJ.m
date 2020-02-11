
function J = leftJ(r)
    tolerance = 1e-12;
    angle = norm(r);
    if angle < tolerance
        % If the angle is small, fall back on the series representation
        N = 10;
        J = eye(3);
        pxn = eye(3);
        px = hat(r);
        for n = 1:N
            pxn = pxn*px/(n + 1);    
            J = J + pxn;
        end
    else
        axis = r/angle;

        cph = (1 - cos(angle))/angle;
        sph = sin(angle)/angle;

        J = sph * eye(3) + (1 - sph) * axis * axis' + cph * hat(axis);
    end       
end
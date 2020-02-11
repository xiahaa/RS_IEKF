
function Jinv = leftJinv(r)
    tolerance = 1e-12;
    phi = r;
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation        
        Jinv = eye(3);
        pxn = eye(3);
        px = hat(phi);
        for n = 1:10
            pxn = pxn * px/n;
            Jinv = Jinv + lut_bernoullinumber(n) * pxn;
        end
    else
        axis = phi/norm(phi);
        ph_2 = 0.5*ph;

        Jinv =   ph_2 * cot(ph_2)* eye(3)...
               + (1 - ph_2 * cot(ph_2))* axis * axis'...
               - ph_2 * hat(axis);
    end   
end
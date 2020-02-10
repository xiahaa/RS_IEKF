function jac = jac_f3(a,b,c,id)
    %% check ok
    if id == 1
        jac = cross(b,c)';
    elseif id == 2
        jac = cross(c,a)';
    elseif id == 3
        jac = cross(a,b)';
    end
end
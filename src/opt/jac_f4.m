function jac = jac_f4(a1,a2,id)
    %% check ok
    if id == 1
        jac = -skewm(a2);
    else
        jac = skewm(a1);
    end
end
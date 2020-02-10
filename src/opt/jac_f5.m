function jac = jac_f5(d)
    %% check ok.
    s1 = norm(d)^3;
    jac = d'*d*eye(3)-d*d';
    jac = jac ./ s1;
end
function Jinv = approxRightJinv(v)
    Jinv = eye(3)+hat(v).*0.5;
end
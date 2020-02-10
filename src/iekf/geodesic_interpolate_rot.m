function varargout = geodesic_interpolate_rot(X, inde, localstamp, t)
% FUNCTION_NAME - geodesic interpolation of rotation matrices.
%
% Syntax:  varargout = geodesic_interpolate_rot(X, inde, localstamp, t)
%
% Inputs:
%    X - staets.
%    inde - index.
%    localstamp - imu timestamp.
%    t - feature timestamp.
%
% Outputs:
%    (1) - interpolated rotation.
%    (2) - relative time.
%    (3) - index of R0.
%    (4) - index of R1.
%    (5) - R0.
%    (6) - R1.
%    (7) - logSO3(R0'*R1)*dt;.
%    (8) - 1/(t_1-t_0).
%
% Author: Xiao Hu
% Technical University of Denmark
% email: xiahaa@space.dtu.dk
% Jan 2020; Last revision: 31-Jan-2020
%------------- BEGIN CODE --------------
    j = find(localstamp>t,1);
    i = j - 1;
    if i < 0, i = 1; warning('extra-interpolation!!! not accurate!!!'); end
    dt = (t - localstamp(i))/(localstamp(j)-localstamp(i));
    s = 1 / (localstamp(j)-localstamp(i));
    i = length(localstamp)+1 - i;
    j = length(localstamp)+1 - j;
    R0 = reshape(X(inde.group(i*9-8:i*9)),3,3);
    R1 = reshape(X(inde.group(j*9-8:j*9)),3,3);
    xit = logSO3(R0'*R1)*dt;
    Rt = R0*expSO3(xit);
    varargout{1} = Rt;
    varargout{2} = dt;
    varargout{3} = i;
    varargout{4} = j;
    varargout{5} = R0;
    varargout{6} = R1;
    varargout{7} = xit;
    varargout{8} = s;
end
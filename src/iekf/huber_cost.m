function [cost, costd] = huber_cost(yhat, delta)
%FUNCTION_NAME - huber_cost.
%
% Syntax:  [cost, costd] = huber_cost(yhat, delta)
%
% Inputs:
%    yhat - error.
%    delta - threshold.
%
% Outputs:
%    cost - robust cost.
%    costd - derivative of robust cost.
%
%
% Author: Xiao Hu
% Technical University of Denmark
% email: xiahaa@space.dtu.dk
% Jan 2020; Last revision: 31-Jan-2020
%------------- BEGIN CODE --------------
    id1 = abs(yhat) <= delta;
    cost = yhat;
    cost(id1) = 0.5*(yhat(id1)).^2;
    id2 = yhat > delta;
    id3 = yhat < -delta;
    cost(id2 | id3) = delta * abs(yhat(id2 | id3)) - 0.5*delta*delta;
    costd = cost;
    costd(id1) = yhat(id1);
    costd(id2) = delta;
    costd(id3) = -delta;
end
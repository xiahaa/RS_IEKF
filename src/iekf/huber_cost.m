function [cost, costd] = huber_cost(yhat, delta)
    id1 = abs(yhat) <= delta;
    cost = yhat;
    cost(id1) = 0.5*(-yhat(id1)).^2;
    id2 = yhat > delta;
    id3 = yhat < -delta;
    cost(id2 | id3) = delta * abs(yhat(id2 | id3)) - 0.5*delta*delta;
    costd = cost;
    costd(id1) = yhat(id1);
    costd(id2) = -delta;
    costd(id3) = delta;
end